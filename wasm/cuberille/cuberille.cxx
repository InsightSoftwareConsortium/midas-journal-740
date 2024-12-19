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

template <typename TImage, typename TMesh, typename TInterpolator>
int cuberilleWithInterpolator(itk::wasm::Pipeline & pipeline, const TImage * image)
{
  using ImageType = TImage;
  using MeshType = TMesh;
  using InterpolatorType = TInterpolator;

  pipeline.get_option("image")->required()->type_name("INPUT_IMAGE");

  double isoSurfaceValue = 1.0;
  pipeline.add_option("--iso-surface-value", isoSurfaceValue, "Value of the iso-surface for which to generate the mesh. Pixels equal to or greater than this value are considered to lie on the surface or inside the resultant mesh.");

  bool quadrilateralFaces = false;
  pipeline.add_flag("--quadrilateral-faces", quadrilateralFaces, "Generate quadrilateral faces instead of triangle faces.");

  bool verticesNotProjectedToIsoSurface = false;
  pipeline.add_flag("--no-projection", verticesNotProjectedToIsoSurface, "Do not project the vertices to the iso-surface.");

  bool imagePixelToCellData = false;
  pipeline.add_flag("--image-pixel-to-cell-data", imagePixelToCellData, "Whether the adjacent input pixel value should be saved as cell data in the output mesh.");

  double projectVertexSurfaceDistanceThreshold = 0.01;
  pipeline.add_option("--surface-distance-threshold", projectVertexSurfaceDistanceThreshold, "Threshold for the distance from the rface during vertex projection in pixel units. Smaller is smoother but takes longer.");

  double projectVertexStepLength = -1.0;
  pipeline.add_option("--step-length", projectVertexStepLength, "Initial step length for vertex projection in physical units. Default is max spacing * 0.35.");

  double projectVertexStepLengthRelaxationFactor = 0.75;
  pipeline.add_option("--step-relaxation-factor", projectVertexStepLengthRelaxationFactor, "The step length relaxation factor during vertex projection. The step length is multiplied by this factor each iteration to allow convergence, [0.0, 1.0].")->check(CLI::Range(0.0, 1.0));

  uint32_t projectVertexMaximumNumberOfSteps = 150;
  pipeline.add_option("--maximum-steps", projectVertexMaximumNumberOfSteps, "The maximum number of steps used during vertex projection.");

  itk::wasm::OutputMesh<MeshType> mesh;
  pipeline.add_option("mesh", mesh, "Output mesh.")->required()->type_name("OUTPUT_MESH");

  ITK_WASM_PARSE(pipeline);

  using InterpolatorBaseType = itk::InterpolateImageFunction<ImageType>;
  using CuberilleType = itk::CuberilleImageToMeshFilter<ImageType, MeshType, InterpolatorBaseType>;

  const auto cuberille = CuberilleType::New();

  const auto interpolator = InterpolatorType::New();
  interpolator->SetInputImage(image);
  cuberille->SetInterpolator(interpolator);

  cuberille->SetInput(image);
  cuberille->SetIsoSurfaceValue(isoSurfaceValue);
  cuberille->SetGenerateTriangleFaces(!quadrilateralFaces);
  cuberille->SetProjectVerticesToIsoSurface(!verticesNotProjectedToIsoSurface);
  cuberille->SetSavePixelAsCellData(imagePixelToCellData);
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

  std::string interpolationMethod = "linear";
  pipeline.add_option("--interpolator", interpolationMethod, "Interpolation method to use for the image. Valid values: linear, bspline, windowed-sync.")->check(CLI::IsMember({"linear", "bspline", "windowed-sinc"}));

  ITK_WASM_PRE_PARSE(pipeline);

  if (interpolationMethod == "linear")
  {
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType>;
    return cuberilleWithInterpolator<ImageType, MeshType, InterpolatorType>(pipeline, image);
  }
  else if (interpolationMethod == "bspline")
  {
    using InterpolatorType = itk::BSplineInterpolateImageFunction<ImageType>;
    return cuberilleWithInterpolator<ImageType, MeshType, InterpolatorType>(pipeline, image);
  }
  else if (interpolationMethod == "windowed-sinc")
  {
    constexpr unsigned int VRadius = 3;
    using InterpolatorType = itk::WindowedSincInterpolateImageFunction<ImageType, VRadius>;
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
  itk::wasm::Pipeline pipeline("cuberille", "Create a mesh from an image via cuberille implicit surface polygonization.", argc, argv);

  return itk::wasm::SupportInputImageTypes<PipelineFunctor,
    uint8_t, int8_t, uint16_t, int16_t, uint32_t, int32_t, float, double>
    ::Dimensions<2U, 3U>("image", pipeline);
}
