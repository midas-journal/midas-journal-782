#include <iostream>
#include <stdlib.h>

#include "itkImageFileReader.h"
#include "itkVectorCastImageFilter.h"
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
  typedef itk::RGBPixel <unsigned char >           RGBPixelType;
  typedef itk::Image< RGBPixelType, Dimension >    RGBImageType;
  typedef itk::Vector< float, Channels >           VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension > VectorImageType;

  // Reader
  typedef itk::ImageFileReader< RGBImageType > FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName( input_file );

  // Cast Filter
  typedef itk::VectorCastImageFilter
    < RGBImageType, VectorImageType > CastFilterType;
  CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput( reader->GetOutput() );

  // Anisotropic Diffusion Filter
  typedef itk::VectorGradientAnisotropicDiffusionImageFilter
    < VectorImageType, VectorImageType >  DiffusionFilterType;
  DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
  diffusion->SetNumberOfIterations( iterations );
  diffusion->SetConductanceParameter( conductance );
  diffusion->SetTimeStep( time_step );
//  diffusion->SetUseImageSpacing( true );
  diffusion->SetInput( caster->GetOutput() );

  // Writer
  typedef itk::ImageFileWriter< VectorImageType > FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( output_file );
  writer->SetInput( diffusion->GetOutput() );

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
