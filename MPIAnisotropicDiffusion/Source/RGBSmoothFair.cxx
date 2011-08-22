#include <iostream>
#include <stdlib.h>

#include "itkImageFileReader.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkImageFileWriter.h"


int main( int argc, char *argv[] )
{
  // Check Command Line
  if (argc < 6 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " input_file conductance iterations time_step " << std::endl;
    std::cerr << " output_file" << std::endl;
    return EXIT_FAILURE;
    }

  // Command Line Parameters
  const char *       input_file      =       argv[1];
  const double       conductance     = atof( argv[2] );
  const unsigned int iterations      = atoi( argv[3] );
  const double       time_step       = atof( argv[4] );
  const char *       output_file     =       argv[5];

  // Image Types
  const unsigned int Dimension = 3;
  const unsigned int Channels = 3;
  typedef itk::Vector< float, Channels >           VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension > VectorImageType;
  typedef VectorImageType::RegionType              RegionType;
  typedef itk::Image< float, Dimension >           ScalarImageType;  
 
  // Reader
  typedef itk::ImageFileReader< VectorImageType > FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName( input_file );
  reader->UpdateOutputInformation();

  VectorImageType::Pointer input_image = reader->GetOutput();
  RegionType  input_region  = input_image->GetLargestPossibleRegion(); 

  // Working Image
  VectorImageType::Pointer working_image = reader->GetOutput(); 

  // Gradient Magnitude Filter
  typedef itk::VectorGradientMagnitudeImageFilter
    < VectorImageType > GradientMagnitudeFilterType;
  GradientMagnitudeFilterType::Pointer 
    gradient = GradientMagnitudeFilterType::New();
  gradient->SetUseImageSpacingOn();
  gradient->SetUsePrincipleComponentsOn();
  gradient->SetInput( working_image );

  unsigned int iteration = 0;
  while ( iteration < iterations )
    {

    try
      {
      gradient->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      }

    // Get the average gradient magnitute.
    typedef itk::ImageRegionConstIterator< ScalarImageType >  IteratorType;
    IteratorType it( gradient->GetOutput(), input_region );
    double sum = 0.0;
    unsigned long int count = input_region.GetNumberOfPixels();
    double avg;
    it.GoToBegin();
    while ( !it.IsAtEnd() )
      {
      sum += it.Get();
      ++it;
      }
     avg = sum / count;
     std::cerr << "Iteration " << iteration
       << ", average gradient magnitude is " << avg << std::endl;

     // Anisotropic Diffusion Filter
     typedef itk::VectorGradientAnisotropicDiffusionImageFilter
       < VectorImageType, VectorImageType >  DiffusionFilterType;
     DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
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
     gradient->SetInput( working_image );
     ++iteration;
    }


  // Writer
  typedef itk::ImageFileWriter< VectorImageType > FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( output_file );
  writer->SetInput( working_image );

  try 
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }


  return EXIT_SUCCESS;
}
