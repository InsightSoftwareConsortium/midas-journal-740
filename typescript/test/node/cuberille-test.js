import test from 'ava'
import path from 'path'

import { readImageNode } from '@itk-wasm/image-io'
import { cuberilleNode } from '../../dist/index-node.js'
import { testInputPath } from './common.js'

test('Test cuberille', async t => {
  const testInputFilePath = path.join(testInputPath, 'fuel.mha')

  const image = await readImageNode(testInputFilePath)
  const { mesh } = await cuberilleNode(image)

  t.assert(mesh, 'cuberille did not return a mesh')
  t.is(mesh.numberOfPoints, 6274)
  t.is(mesh.points.length, 18822)
  t.is(mesh.numberOfCells, 12544)
  t.is(mesh.cells.length, 62720)
})
