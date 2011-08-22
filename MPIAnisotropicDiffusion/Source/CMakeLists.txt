INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

PROJECT( MPIVisibleWoman )

cmake_minimum_required( VERSION 2.8 )

ENABLE_TESTING()
INCLUDE(Dart)

SET(Required_Packages
ITK
)

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

INCLUDE( FindMPI )
INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )

SET( ITK_LIBS ${ITK_LIBS} ITKCommon )
SET( ITK_LIBS ${ITK_LIBS} ITKIO )

SET( MPI_LIBS ${MPI_LIBS} "${MPI_LIBRARY}" )
SET( MPI_LIBS ${MPI_LIBS} "${MPI_EXTRA_LIBRARY}" )


ADD_EXECUTABLE( MPIVectorGradientAnisotropicDiffusionImageFilter
  MPIVectorGradientAnisotropicDiffusionImageFilter.cxx )
TARGET_LINK_LIBRARIES ( MPIVectorGradientAnisotropicDiffusionImageFilter
  ${ITK_LIBS} ${MPI_LIBS} )

ADD_EXECUTABLE(  RGBSmooth RGBSmooth.cxx )
TARGET_LINK_LIBRARIES ( RGBSmooth  ${ITK_LIBS} )

ADD_EXECUTABLE(  RGBSmoothFair RGBSmoothFair.cxx )
TARGET_LINK_LIBRARIES ( RGBSmoothFair  ${ITK_LIBS} )

ADD_TEST(TestSerial ./RGBSmoothFair
  ${CMAKE_SOURCE_DIR}/VisibleWomanEye.mhd 1.0 10 0.015 Serial.mhd )

ADD_TEST(TestMPI ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} 
  ./MPIVectorGradientAnisotropicDiffusionImageFilter ${MPIEXEC_POSTFLAGS} 
  ${CMAKE_SOURCE_DIR}/VisibleWomanEye.mhd MPI.mhd 1.0 10 0.015 )

ADD_TEST(TestCompare ImageCompare MPI.mhd Serial.mhd )

