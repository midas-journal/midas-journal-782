#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "itkVector.h"

#include "itkImageFileReader.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkPasteImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageIORegion.h"

#include "itkImageRegionNonCubeSplitter.h"
#include "itkImageRegionSplitter.h"


#define TAG_SUM     0
#define TAG_AVERAGE 1
#define TAG_BORDER  2
#define TAG_PIECE   3

int main( int argc, char *argv[] )
{

  const char *       in_file_name    =       argv[1];
  const char *       out_file_name   =       argv[2];
  const double       conductance     = atof( argv[3] );
  const unsigned int iterations      = atoi( argv[4] );
  const double       time_step       = atof( argv[5] );

  const unsigned int Dimension = 3;
  const unsigned int Channels  = 3;
  typedef itk::Vector <float, Channels >           VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension > VectorImageType;
  typedef VectorImageType::RegionType              RegionType;
  typedef VectorImageType::IndexType               IndexType;
  typedef VectorImageType::SizeType                SizeType;
  typedef VectorImageType::PointType               PointType;
  typedef VectorImageType::SpacingType             SpacingType;
  typedef itk::Image< float, Dimension >           ScalarImageType;
  
  int mpi_rank;
  int mpi_size;

  // Initialise MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size ); 

  typedef itk::PasteImageFilter < 
    VectorImageType, VectorImageType, VectorImageType > PasteType;

  // Reader
  typedef itk::ImageFileReader< VectorImageType > FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName( in_file_name );
  reader->UpdateOutputInformation();

  VectorImageType::Pointer input_image = reader->GetOutput();
  PointType   input_origin  = input_image->GetOrigin();
  RegionType  input_region  = input_image->GetLargestPossibleRegion();
  IndexType   input_index   = input_region.GetIndex();
  SizeType    input_size    = input_region.GetSize();
  SpacingType input_spacing = input_image->GetSpacing();

  // Region Splitter
  typedef itk::ImageRegionSplitter< Dimension > SplitterType;
  SplitterType::Pointer splitter = SplitterType::New();

  // Get The splits
  std::vector < RegionType > split_regions( mpi_size );
  std::vector < RegionType > padded_regions( mpi_size );
  for ( int split = 0; split < mpi_size; ++split )
    {
    split_regions[ split ] 
      = splitter->GetSplit( split, mpi_size, input_region );
    
    padded_regions[ split ] = split_regions[ split ];
    padded_regions[ split ].PadByRadius( 1 );
    padded_regions[ split ].Crop( input_region );
    
    } // end for split

  // We need to know who our neighbors are 
  // and the regions of overlap. 
  //
  // This will be a "two way street." That is we will send 
  // the part of our split that overlaps our neighbor's 
  // padded region and we will receive the part of our 
  // padded region that overlaps our neighbor's split.

  // Our send regions
  std::vector< RegionType > com_snd_regions( mpi_size );
  for ( int split = 0; split < mpi_size; ++split )
    {
    com_snd_regions[ split ] = split_regions[ mpi_rank ];
    bool good_crop = ( com_snd_regions[ split ].Crop( padded_regions[ split ] ));
    // We don't need to talk to ourself.   
    if (( split == mpi_rank ) || ( !good_crop ))
      {
      SizeType  com_size;
      com_size.Fill( 0 );
      com_snd_regions[ split ].SetSize(  com_size );
      }        
    } // end for split

  // Our receive regions
  std::vector< RegionType > com_rcv_regions( mpi_size );
  for ( int split = 0; split < mpi_size; ++split )
    {
    com_rcv_regions[ split ] = padded_regions[ mpi_rank ];
    bool good_crop = (com_rcv_regions[ split ].Crop( split_regions[ split ] ));
    // We don't need to talk to ourself.   
    if (( split == mpi_rank ) || ( !good_crop ))
      {
      SizeType  com_size;
      com_size.Fill( 0 );
      com_rcv_regions[ split ].SetSize(  com_size );
      }        
    } // end for split
  
  // We need arrays of smart pointers to 
  // images for the send and receive buffers.
  std::vector< VectorImageType::Pointer > send_images( mpi_size );
  std::vector< VectorImageType::Pointer > recv_images( mpi_size );
  
  // We must create the receive images before we can get any data.
  for ( int split = 0; split < mpi_size; ++split )
    {
    if ( com_rcv_regions[ split ].GetNumberOfPixels() != 0 )
      {
      recv_images[ split ] = VectorImageType::New();
      recv_images[ split ]->SetOrigin( input_origin );
      recv_images[ split ]->SetRegions( com_rcv_regions[ split ] );
      recv_images[ split ]->Allocate();
      }
    }

  // Extractor Type
  typedef itk::ExtractImageFilter< 
    VectorImageType, VectorImageType > ExtractorType;

  // Padded Extractor
  ExtractorType::Pointer pad_extractor = ExtractorType::New();
  pad_extractor->SetExtractionRegion( padded_regions[ mpi_rank] );
  pad_extractor->SetInput( reader->GetOutput() );

  try
    {
    pad_extractor->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught updating padded extractor!" << std::endl;
    std::cerr << excep << std::endl;
    }

  // Working Image
  VectorImageType::Pointer working_image = pad_extractor->GetOutput();
  working_image->DisconnectPipeline();
  
   unsigned int iteration = 0;
   while ( iteration < iterations )
     {

     // Gradient Magnitude Filter
     typedef itk::VectorGradientMagnitudeImageFilter
       < VectorImageType > GradientMagnitudeFilterType;
     GradientMagnitudeFilterType::Pointer 
       gradient = GradientMagnitudeFilterType::New();
     gradient->SetUseImageSpacingOn();
     gradient->SetUsePrincipleComponentsOn();
     gradient->SetInput( working_image );
      
     try
       {
       gradient->Update();
       }
     catch( itk::ExceptionObject & excep )
       {
       std::cerr << "Exception caught !" << std::endl;
       std::cerr << excep << std::endl;
       }

     // Get the average gradient magnitute of ALL splits.
     typedef itk::ImageRegionConstIterator< ScalarImageType >  IteratorType;
     IteratorType it( gradient->GetOutput(), split_regions[ mpi_rank ] );
     double sum = 0.0;
     unsigned long int count = split_regions[ mpi_rank ].GetNumberOfPixels();
     double avg;
     it.GoToBegin();
     while ( !it.IsAtEnd() )
       {
       sum += it.Get();
       ++it;
       }
     if ( mpi_rank != 0 )
       {
       MPI_Ssend( &sum, 1, MPI_DOUBLE, 0, TAG_SUM,     MPI_COMM_WORLD );
       MPI_Recv( &avg, 1, MPI_DOUBLE, 0, TAG_AVERAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
       }
     else
       {
       for ( int split = 1; split < mpi_size; ++split )
         {
         double other_sum;
         unsigned long int other_count = split_regions[ split ].GetNumberOfPixels();;
         MPI_Recv( &other_sum, 1, MPI_DOUBLE, split, TAG_SUM, 
           MPI_COMM_WORLD, MPI_STATUS_IGNORE );
         sum   += other_sum;
         count += other_count;
         }
       avg = sum / count;
       std::cerr << "Iteration " << iteration 
         << ", average gradient magnitude is " << avg << std::endl;
       for ( int split = 1; split < mpi_size; ++split )
         {
         MPI_Ssend( &avg, 1, MPI_DOUBLE, split, TAG_AVERAGE, MPI_COMM_WORLD );
         }
       }

     typedef itk::VectorGradientAnisotropicDiffusionImageFilter
        < VectorImageType, VectorImageType > DiffusionType;
      DiffusionType::Pointer diffusion = DiffusionType::New();
      diffusion->SetNumberOfIterations( 1 );
      diffusion->SetConductanceParameter( conductance );
      diffusion->SetTimeStep( time_step );
      diffusion->SetFixedAverageGradientMagnitude( avg );
      diffusion->SetInput( working_image );

      try
        {
        diffusion->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        }

      // Disconnect the working image from the pipeline.
      working_image = diffusion->GetOutput();
      working_image->DisconnectPipeline();

      // Prepare the updated boundaries
      for ( int split = 0; split < mpi_size; ++split )
        {
        if ( com_snd_regions[ split ].GetNumberOfPixels() != 0 )
          {
          // Send Extractor
          ExtractorType::Pointer send_extractor = ExtractorType::New();
          send_extractor->SetInput( working_image );
          send_extractor->SetExtractionRegion( com_snd_regions[ split ] );
          try
            {
            send_extractor->Update();
            }
          catch( itk::ExceptionObject & excep )
            {
            std::cerr << "Exception caught while updating Send Extractor!" << std::endl;
            std::cerr << excep << std::endl;
            }
          send_images[ split ] = send_extractor->GetOutput();
          send_images[ split ]->DisconnectPipeline();
          } // end prepare send neighbors
        } // end prepare send for split

      // MPI Send/Receive
      for ( int split = 0; split < mpi_size; ++split )
        {
        if ( mpi_rank < split )
          {
          if ( com_snd_regions[ split ].GetNumberOfPixels() != 0 )
            {
            MPI_Ssend( send_images[ split ]->GetBufferPointer(), 
              com_snd_regions[ split ].GetNumberOfPixels() * Channels,
              MPI_FLOAT, split, TAG_BORDER, MPI_COMM_WORLD );
            }
          if ( com_rcv_regions[ split ].GetNumberOfPixels() != 0 )
            {
            MPI_Recv( recv_images[ split ]->GetBufferPointer(),
              com_rcv_regions[ split ].GetNumberOfPixels() * Channels,
              MPI_FLOAT, split, TAG_BORDER, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            }
          }
        else
          {
          if ( com_rcv_regions[ split ].GetNumberOfPixels() != 0 )
            {
            MPI_Recv( recv_images[ split ]->GetBufferPointer(),
              com_rcv_regions[ split ].GetNumberOfPixels() * Channels,
              MPI_FLOAT, split, TAG_BORDER, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            }
          if ( com_snd_regions[ split ].GetNumberOfPixels() != 0 )
            {
            MPI_Ssend( send_images[ split ]->GetBufferPointer(), 
              com_snd_regions[ split ].GetNumberOfPixels() * Channels,
              MPI_FLOAT, split, TAG_BORDER, MPI_COMM_WORLD );
            }
          }
        } // end send/receive for split

      // Use the updated boundaries
      for ( int split = 0; split < mpi_size; ++split )
        {
        if ( com_rcv_regions[ split ].GetNumberOfPixels() != 0 )
          {
          // Paste the receive image into the working image.
          PasteType::Pointer paste = PasteType::New();
          paste->SetDestinationImage( working_image );
          paste->SetSourceImage( recv_images[ split ] );
          paste->SetDestinationIndex( com_rcv_regions[ split ].GetIndex() );
          paste->SetSourceRegion( com_rcv_regions[ split ] );
          try
            {
            paste->Update();
            }
          catch( itk::ExceptionObject & excep )
            {
            std::cerr << "Exception caught while updating receive paste!" << std::endl;
            std::cerr << excep << std::endl;
            }
          working_image = paste->GetOutput();
          working_image->DisconnectPipeline();

          } // end use neighbors
        } // end use for split

      ++iteration;

     } // end of iteration loop

  /////////////////////////////////
  // Start the final output here //
  /////////////////////////////////

   ExtractorType::Pointer final_extractor = ExtractorType::New();
   final_extractor->SetInput( working_image );
   final_extractor->SetExtractionRegion( split_regions[ mpi_rank ] );

   try
     {
     final_extractor->Update();
     }
   catch( itk::ExceptionObject & excep )
     {
     std::cerr << "Exception caught while updating final extractor!" << std::endl;
     std::cerr << excep << std::endl;
     }

  // Each process, other than 0, sends it's image to process 0.
  if ( mpi_rank != 0 )
    {
    MPI_Ssend( 
      final_extractor->GetOutput()->GetBufferPointer(),
      split_regions[ mpi_rank ].GetNumberOfPixels() * Channels,
      MPI_FLOAT, 0, TAG_PIECE, MPI_COMM_WORLD );
    }
  else
    {

    
    for ( int split = 0; split < mpi_size; ++split )
      {

      
      // Get the receive image region.
      RegionType recv_region = split_regions[ split ];
        //= splitter->GetSplit( split, mpi_size, input_region );
      IndexType write_index = recv_region.GetIndex();
      SizeType  write_size  = recv_region.GetSize();

      // Make the receive image.
      VectorImageType::Pointer recv_image = VectorImageType::New();
      recv_image->SetOrigin( input_origin );
      recv_image->SetRegions( recv_region );
      if ( split != 0 )
        {
        recv_image->Allocate();

        // Receive the image.
        MPI_Recv( recv_image->GetBufferPointer(), 
          recv_region.GetNumberOfPixels() * Channels, 
          MPI_FLOAT, split, TAG_PIECE, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }

      // Paste
      PasteType::Pointer paste = PasteType::New();
      paste->SetDestinationImage( reader->GetOutput() );
      if ( split == 0 )
        {
        paste->SetSourceImage( final_extractor->GetOutput() );
        }
      else
        {
        paste->SetSourceImage( recv_image );
        }
      paste->SetDestinationIndex( write_index );
      paste->SetSourceRegion( recv_region );
        
      // Writer
      typedef itk::ImageFileWriter< VectorImageType > FileWriterType;
      FileWriterType::Pointer writer = FileWriterType::New();
      writer->SetFileName( out_file_name );
      writer->SetInput( paste->GetOutput() );

      itk::ImageIORegion write_region( Dimension );
      for ( unsigned int dim = 0; dim < Dimension; ++dim )
        {
        write_region.SetIndex( dim, write_index[ dim ] );
        write_region.SetSize(  dim, write_size[  dim ] );
        }
      writer->SetIORegion( write_region );

      try
        {
         std::cerr << "Piece " << split << 
           " updating writer." << std::endl;
        writer->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        }
      } // end split loop
    } // end rank 0
   
  MPI_Finalize();


  return EXIT_SUCCESS;
}
