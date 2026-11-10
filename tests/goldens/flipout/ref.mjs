import fs from 'node:fs'
import factory from './nxr_compute.js'
const mod = await factory()
const M = JSON.parse(fs.readFileSync('meshes.json'))
const out = {}
for (const [name, m] of Object.entries(M)) {
  const ctx = new mod.Manifold(new Float64Array(m.vertices.flat()), new Int32Array(m.faces.flat()))
  const paths = []
  for (const [a, b] of m.pairs) {
    if (a === b) continue
    const r = ctx.tracePath(a, b)
    paths.push({ from: a, to: b, length: r.length, exact: ctx.exactDistance(a, b), nPoints: r.nPoints, positions: Array.from(r.positions), vertices: Array.from(r.vertices), weights: Array.from(r.weights) })
  }
  out[name] = { ...m, paths }
  delete out[name].pairs
  if (name.startsWith('fsaverage')) { delete out[name].vertices; delete out[name].faces }
}
fs.writeFileSync('flipout-reference.json', JSON.stringify(out))
console.log(Object.fromEntries(Object.entries(out).map(([k, v]) => [k, v.paths.length])))
