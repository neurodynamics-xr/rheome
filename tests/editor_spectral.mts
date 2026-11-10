/**
 * The editor's own CPU functions on one chunk, for tSpectralRoot's parity: the load (chunkFft.ts: less the row mean, mirror-padded),
 * fftCpu (fft.oracle), N·mean back in bin 0, the one-sided spectrum held as float32 (as on the device), bandGains, spectralColumnCpu
 * and spectralAnalyticCpu at the asked columns.  npx tsx editor_spectral.mts <editor> <dir>  (editor: the Cortical Flow checkout; dir: rows.f32 [C × w] row-major, meta.json)
 */
import { readFileSync, writeFileSync } from 'node:fs'
import { join } from 'node:path'

const R = process.argv[2]!, E = `${R}/editor/src`
const { fftCpu, nextPow2 } = await import(`${E}/engine/toolbox/primitives/fft/fft.oracle.ts`)
const { meyerBand } = await import(`${E}/engine/toolbox/wavelets/dyadicFrame.ts`)
/* modeView.ts imports the device (three/tsl) at load: its three CPU functions are taken from its SOURCE, verbatim, and compiled alone */
const { transformSync } = await import(`${R}/node_modules/esbuild/lib/main.js`)
const src = readFileSync(`${E}/editor/engine/gpu/modeView.ts`, 'utf8')
const fn = (name: string) => { const a = src.indexOf(`export function ${name}(`); const b = src.indexOf('\n}\n', a); if (a < 0 || b < 0) throw new Error(name); return src.slice(a, b + 3) }
;(globalThis as any).__meyerBand = meyerBand
const js = transformSync(['spectralColumnCpu', 'bandGains', 'spectralAnalyticCpu'].map(fn).join('\n'), { loader: 'ts', format: 'esm' }).code
const { bandGains, spectralColumnCpu, spectralAnalyticCpu } = await import('data:text/javascript;base64,' + Buffer.from('const meyerBand = globalThis.__meyerBand;\n' + js).toString('base64'))

const dir = process.argv[3]!
const meta = JSON.parse(readFileSync(join(dir, 'meta.json'), 'utf8')) as { C: number; w: number; halo: number; rate: number; band: [number, number]; cols: number[] }
const buf = readFileSync(join(dir, 'rows.f32')), rows = new Float32Array(buf.buffer, buf.byteOffset, buf.byteLength / 4)
const { C, w, halo } = meta, N = nextPow2(w), bins = N / 2 + 1
const X = new Float32Array(C * bins * 2)
for (let c = 0; c < C; c++) {
  let m = 0
  for (let j = 0; j < w; j++) m += rows[c * w + j]!
  m /= w
  const re = new Float64Array(N), im = new Float64Array(N)
  for (let i = 0; i < N; i++) { const src = i < w ? i : Math.max(0, 2 * (w - 1) - i); re[i] = rows[c * w + src]! - m }
  const s = fftCpu(re, im)
  for (let k = 0; k < bins; k++) { X[2 * (c * bins + k)] = s.re[k]! + (k === 0 ? N * m : 0); X[2 * (c * bins + k) + 1] = s.im[k]! }
}
const H = bandGains(bins, N, meta.rate, meta.band)
const out: number[] = []
for (const col of meta.cols) {
  const n = halo + col
  const plain = spectralColumnCpu(X, C, N, n, null), { re, im } = spectralAnalyticCpu(X, C, N, n, H)
  for (let c = 0; c < C; c++) out.push(plain[c]!, re[c]!, im[c]!)
}
writeFileSync(join(dir, 'cols.f64'), Buffer.from(new Float64Array(out).buffer))
writeFileSync(join(dir, 'H.f64'), Buffer.from(Float64Array.from(H).buffer))
writeFileSync(join(dir, 'X.f32'), Buffer.from(X.buffer))
console.log(JSON.stringify({ N, bins, cols: meta.cols.length }))
