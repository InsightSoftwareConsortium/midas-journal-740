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

#include "itkPipeline.h"
#include "itkInputImage.h"
#include "itkOutputMesh.h"
#include "itkSupportInputImageTypes.h"

template <typename TImage>
int cuberille(itk::wasm::Pipeline &pipeline, const TImage * image)
{
  using ImageType = TImage;
  constexpr unsigned int Dimension = ImageType::ImageDimension;
  using PixelType = typename ImageType::PixelType;  using MeshType = itk::Mesh<PixelType, Dimension>;


  pipeline.get_option("image")->required()->type_name("INPUT_IMAGE");

  bool quadrilateralFaces = false;
  pipeline.add_flag("--quadrilateral-faces", quadrilateralFaces, "Generate quadrilateral faces instead of triangle faces.");

  bool verticesNotProjectedToIsoSurface = false;
  pipeline.add_flag("--vertices-not-projected-to-iso-surface", verticesNotProjectedToIsoSurface, "Do not project the vertices to the iso-surface.");

  bool imagePixelToCellData = false;
  pipeline.add_flag("--image-pixel-to-cell-data", imagePixelToCellData, "Whether the adjacent input pixel value should be saved as cell data in the output mesh.");

  double projectVertexSurfaceDistanceThreshold = 0.5;
  pipeline.add_option("--project-vertex-surface-distance-threshold", projectVertexSurfaceDistanceThreshold, "Threshold for the distance from the iso-surface during vertex projection in pixel units. Smaller is smoother but takes longer.");

  double projectVertexStepLength = -1.0;
  pipeline.add_option("--project-vertex-step-length", projectVertexStepLength, "Initial step length for vertex projection in physical units. Default is max spacing * 0.25.");

  double projectVertexStepLengthRelaxationFactor = 0.95;
  pipeline.add_option("--project-vertex-step-length-relaxation-factor", projectVertexStepLengthRelaxationFactor, "The step length relaxation factor during vertex projection. The step length is multiplied by this factor each iteration to allow convergence, [0.0, 1.0].");

  uint32_t projectVertexMaximumNumberOfSteps = 50;
  pipeline.add_option("--project-vertex-maximum-number-of-steps", projectVertexMaximumNumberOfSteps, "The maximum number of steps used during vertex projection.");

  bool generateQuadrilateralFaces = false;
  pipeline.add_flag("--generate-quadrilateral-faces", generateQuadrilateralFaces, "Generate quadrilateral faces instead of triangle faces.");

  itk::wasm::OutputMesh<MeshType> mesh;
  pipeline.add_option("mesh", mesh, "Output mesh.")->required()->type_name("OUTPUT_MESH");

  ITK_WASM_PARSE(pipeline);

  // Pipeline code goes here

  return EXIT_SUCCESS;
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
