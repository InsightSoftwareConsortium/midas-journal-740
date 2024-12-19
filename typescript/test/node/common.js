import path from 'path'
import fs from 'fs'

export const testInputPath = path.resolve('..', 'test', 'Input')
export const testBaselinePath = path.resolve('..', 'test', 'Baseline')
export const testOutputPath = path.resolve('..', 'test', 'Output', 'typescript')
fs.mkdirSync(testOutputPath, { recursive: true })