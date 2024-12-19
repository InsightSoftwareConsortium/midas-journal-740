/*=========================================================================

 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkMesh.h"

#include "itkMeshJSON.h"
#include "itkPipeline.h"
#include "itkInputImage.h"
#include "itkOutputMesh.h"
#include "itkSupportInputImageTypes.h"

#include "itkCuberilleImageToMeshFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkFlatStructuringElement.h"

template <typename TImage, typename TMesh, typename TInterpolator>
int cuberilleWithInterpolator(itk::wasm::Pipeline & pipeline, const TImage * image)
{
  using ImageType = TImage;
  using MeshType = TMesh;
  using InterpolatorType = TInterpolator;

  pipeline.get_option("image")->required()->type_name("INPUT_IMAGE");

  double maximumRMSError = 0.07;
  pipeline.add_flag("--maximum-rms-error", maximumRMSError, "Maximum root mean square error for the anti-aliasing filter. Smaller is smoother but takes longer.");

  uint32_t numberOfIterations = 1000;
  pipeline.add_flag("--iterations", numberOfIterations, "Number of iterations for the anti-aliasing filter. Larger is smoother but takes longer.");

  double smoothingSigma = 0.95;
  pipeline.add_option("--smoothing-sigma", smoothingSigma, "Sigma for the Gaussian smoothing filter in pixel units.");

  bool noFillHoles = false;
  pipeline.add_flag("--no-fill-holes", noFillHoles, "Do not fill holes in the binary image before anti-aliasing.");

  bool noClosing = false;
  pipeline.add_flag("--no-closing", noClosing, "Do not close the smoothed, anti-aliased image before mesh generation.");

  double isoSurfaceValue = 0.0;
  pipeline.add_flag("--iso-surface-value", isoSurfaceValue, "Value of the iso-surface on the normalized, anti-aliased image for which to generate the mesh. Pixels equal to or greater than this value are considered to lie on the surface or inside the resultant mesh. [-4.0, 4.0]");

  bool quadrilateralFaces = false;
  pipeline.add_flag("--quadrilateral-faces", quadrilateralFaces, "Generate quadrilateral faces instead of triangle faces.");

  double projectVertexSurfaceDistanceThreshold = 0.01;
  pipeline.add_option("--surface-distance-threshold", projectVertexSurfaceDistanceThreshold, "Threshold for the distance from the iso-surface during vertex projection in pixel units. Smaller is smoother but takes longer.");

  double projectVertexStepLength = -1.0;
  pipeline.add_option("--step-length", projectVertexStepLength, "Initial step length for vertex projection in physical units. Default is max spacing * 0.35.");

  double projectVertexStepLengthRelaxationFactor = 0.75;
  pipeline.add_option("--step-relaxation-factor", projectVertexStepLengthRelaxationFactor, "The step length relaxation factor during vertex projection. The step length is multiplied by this factor each iteration to allow convergence, [0.0, 1.0].")->check(CLI::Range(0.0, 1.0));

  uint32_t projectVertexMaximumNumberOfSteps = 150;
  pipeline.add_option("--maximum-steps", projectVertexMaximumNumberOfSteps, "The maximum number of steps used during vertex projection.");

  itk::wasm::OutputMesh<MeshType> mesh;
  pipeline.add_option("mesh", mesh, "Output mesh.")->required()->type_name("OUTPUT_MESH");

  ITK_WASM_PARSE(pipeline);

  using FillholeFilterType = itk::BinaryFillholeImageFilter<ImageType>;
  const auto fillholeFilter = FillholeFilterType::New();
  fillholeFilter->SetInput(image);

  using FloatImageType = itk::Image<float, ImageType::ImageDimension>;
  using AntiAliasFilterType = itk::AntiAliasBinaryImageFilter<ImageType, FloatImageType>;
  const auto antiAliasFilter = AntiAliasFilterType::New();
  if (noFillHoles)
  {
    antiAliasFilter->SetInput(image);
  }
  else
  {
    antiAliasFilter->SetInput(fillholeFilter->GetOutput());
  }
  antiAliasFilter->SetMaximumRMSError(maximumRMSError);
  antiAliasFilter->SetNumberOfIterations(numberOfIterations);

  using SmoothingFilterType = itk::SmoothingRecursiveGaussianImageFilter<FloatImageType, FloatImageType>;
  const auto smoothingFilter = SmoothingFilterType::New();
  smoothingFilter->SetInput(antiAliasFilter->GetOutput());
  typename SmoothingFilterType::SigmaArrayType sigma;
  for (unsigned int i = 0; i < ImageType::ImageDimension; ++i)
  {
    sigma[i] = image->GetSpacing()[i] * smoothingSigma;
  }
  smoothingFilter->SetSigmaArray(sigma);

  using StructuringElementType = itk::FlatStructuringElement<ImageType::ImageDimension>;
  using ClosingFilterType = itk::GrayscaleMorphologicalClosingImageFilter<FloatImageType, FloatImageType, StructuringElementType>;
  const auto closingFilter = ClosingFilterType::New();
  closingFilter->SetInput(smoothingFilter->GetOutput());
  typename StructuringElementType::RadiusType closingRadius;
  closingRadius.Fill(1);
  const auto closingStructuringElement = StructuringElementType::Ball(closingRadius, false);
  closingFilter->SetKernel(closingStructuringElement);

  using InterpolatorBaseType = itk::InterpolateImageFunction<FloatImageType>;
  using CuberilleType = itk::CuberilleImageToMeshFilter<FloatImageType, MeshType, InterpolatorBaseType>;

  const auto cuberille = CuberilleType::New();

  const auto interpolator = InterpolatorType::New();
  cuberille->SetInterpolator(interpolator);

  if (noClosing)
  {
    interpolator->SetInputImage(smoothingFilter->GetOutput());
    cuberille->SetInput(smoothingFilter->GetOutput());
  }
  else
  {
    interpolator->SetInputImage(closingFilter->GetOutput());
    cuberille->SetInput(closingFilter->GetOutput());
  }
  cuberille->SetIsoSurfaceValue(isoSurfaceValue);
  cuberille->SetGenerateTriangleFaces(!quadrilateralFaces);
  cuberille->SetProjectVertexSurfaceDistanceThreshold(projectVertexSurfaceDistanceThreshold);
  if (projectVertexStepLength > 0.0)
  {
    cuberille->SetProjectVertexStepLength(projectVertexStepLength);
  }
  else
  {
    cuberille->SetProjectVertexStepLength(image->GetSpacing().GetVnlVector().max_value() * 0.35);
  }
  cuberille->SetProjectVertexStepLengthRelaxationFactor(projectVertexStepLengthRelaxationFactor);
  cuberille->SetProjectVertexMaximumNumberOfSteps(projectVertexMaximumNumberOfSteps);

  ITK_WASM_CATCH_EXCEPTION(pipeline, cuberille->Update());

  const typename MeshType::ConstPointer outputMesh = cuberille->GetOutput();
  mesh.Set(outputMesh);

  return EXIT_SUCCESS;
}

template <typename TImage>
int cuberille(itk::wasm::Pipeline &pipeline, const TImage * image)
{
  using ImageType = TImage;
  constexpr unsigned int Dimension = ImageType::ImageDimension;
  using PixelType = typename ImageType::PixelType;
  using MeshType = itk::Mesh<PixelType, Dimension>;

  using FloatImageType = itk::Image<float, ImageType::ImageDimension>;

  std::string interpolationMethod = "linear";
  pipeline.add_option("--interpolator", interpolationMethod, "Interpolation method to use for the image. Valid values: linear, bspline, windowed-sync.")->check(CLI::IsMember({"linear", "bspline", "windowed-sinc"}));

  ITK_WASM_PRE_PARSE(pipeline);

  if (interpolationMethod == "linear")
  {
    using InterpolatorType = itk::LinearInterpolateImageFunction<FloatImageType>;
    return cuberilleWithInterpolator<ImageType, MeshType, InterpolatorType>(pipeline, image);
  }
  else if (interpolationMethod == "bspline")
  {
    using InterpolatorType = itk::BSplineInterpolateImageFunction<FloatImageType>;
    return cuberilleWithInterpolator<ImageType, MeshType, InterpolatorType>(pipeline, image);
  }
  else if (interpolationMethod == "windowed-sinc")
  {
    constexpr unsigned int VRadius = 3;
    using InterpolatorType = itk::WindowedSincInterpolateImageFunction<FloatImageType, VRadius>;
    return cuberilleWithInterpolator<ImageType, MeshType, InterpolatorType>(pipeline, image);
  }
  else
  {
    std::ostringstream ostrm;
    ostrm << "Unknown interpolation method: " << interpolationMethod;
    CLI::Error err("Runtime error", ostrm.str(), 1);
    return pipeline.exit(err);
  }
}

template <typename TImage>
class PipelineFunctor
{
public:
  int operator()(itk::wasm::Pipeline &pipeline)
  {
    using ImageType = TImage;

    itk::wasm::InputImage<ImageType> image;
    pipeline.add_option("image", image, "Input image")->type_name("INPUT_IMAGE");

    ITK_WASM_PRE_PARSE(pipeline);

    typename ImageType::ConstPointer imageRef = image.Get();
    return cuberille<ImageType>(pipeline, imageRef);
  }
};

int main(int argc, char * argv[])
{
  itk::wasm::Pipeline pipeline("cuberille", "Anti-alias a label image and create a mesh via cuberille implicit surface polygonization.", argc, argv);

  return itk::wasm::SupportInputImageTypes<PipelineFunctor,
    uint8_t, int8_t, uint16_t, int16_t, uint32_t, int32_t>
    ::Dimensions<2U, 3U>("image", pipeline);
}
