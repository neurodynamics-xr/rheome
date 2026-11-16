/**
 * Put the Studio in a page and show it a SCENE (studio/src/studio/show/scene.ts; spec docs/superpowers/specs/2026-11-10-studio-commands.md).
 *
 *   const studio = await mount(element, { studio: '/', theme: 'light' })
 *   await studio.show(scene)                       // one transaction, one commit: the first frame is the scene; again: reconciled by name
 *   const { measurements } = await studio.measure({ path_mm: { length: 'path' } })
 *   const { png } = await studio.export({ png: { scale: 2, background: 'white' } })
 *
 * Only messages cross the frame, so any origin works. The first show opens the Studio on the scene's loaded surface — or, for a scene with no surface
 * (`stage: false`, its fields on three lines), the scopes window: the scope alone.
 */
/** a scene's bundles and values files named relative to THIS page, made absolute (the Studio resolves nothing against a page it cannot see) */
const absolute = (scene) => {
  const s = structuredClone(scene)
  for (const m of Object.values(s.manifolds ?? {})) if (m?.load?.bundle) {
    const u = new URL(m.load.bundle, location.href); m.load.bundle = u.pathname.endsWith('.json') ? u.href : new URL('manifest.json', u.href.endsWith('/') ? u.href : u.href + '/').href
  }
  for (const f of Object.values(s.fields ?? {})) if (f?.load?.values) f.load.values = new URL(f.load.values, location.href).href
  for (const m of Object.values(s.manifolds ?? {})) if (m?.load?.root) m.load.root = new URL(m.load.root.endsWith('/') ? m.load.root : m.load.root + '/', location.href).href.replace(/\/$/, '')
  for (const f of Object.values(s.fields ?? {})) if (f?.load?.sensors?.root) f.load.sensors.root = new URL(f.load.sensors.root.endsWith('/') ? f.load.sensors.root : f.load.sensors.root + '/', location.href).href.replace(/\/$/, '')
  return s
}

export async function mount(element, { studio = '/', theme = 'light' } = {}) {
  const frame = document.createElement('iframe')
  frame.src = new URL(`index.html?embed=1&theme=${theme}`, new URL(studio, location.href)).href
  frame.allow = 'fullscreen'
  Object.assign(frame.style, { width: '100%', height: '100%', border: '0', display: 'block' })
  element.appendChild(frame)
  const origin = new URL(frame.src).origin, waiting = new Map()
  let n = 0, hello, failed = null
  const helloP = new Promise((ok) => { hello = ok })
  window.addEventListener('message', (e) => {
    if (e.source !== frame.contentWindow || e.origin !== origin) return
    const d = e.data
    if (!d || (d.type !== 'nxr-studio' && d.type !== 'nxr-editor')) return   // nxr-editor: the Studio's builds before its rename (2026-11-16)
    if (d.state === 'hello') hello()
    else if (d.state === 'result' && waiting.has(d.id)) { waiting.get(d.id)(d); waiting.delete(d.id) }
    else if (d.state === 'error') { failed = d.code === 'no-webgpu' ? 'this browser has no WebGPU (use Chrome or Edge)' : d.message; for (const [id, ok] of waiting) { ok({ ok: false, error: { section: 'scene', name: '', message: failed } }); waiting.delete(id) } }
  })
  await helloP
  const ask = (op, body) => new Promise((ok) => {
    if (failed) return ok({ ok: false, error: { section: 'scene', name: '', message: failed } })
    const id = `${op}-${++n}`
    waiting.set(id, ok)
    frame.contentWindow.postMessage({ type: 'nxr-host', op, id, ...body }, origin)
  })
  return {
    frame,
    /** show a scene: { ok } or { ok: false, error: { section, name, message } } */
    show: (scene) => ask('show', { scene: absolute(scene) }),
    /** measure what the shown scene made: { ok, measurements: { name: number } } */
    measure: (request) => ask('measure', { request }),
    /** export the shown scene: { ok, png (base64), svg, meta } — the Studio's figure export of the Stage */
    export: (request) => ask('export', { request }),
    /** a cortex panel of the shown bundle (the request: the viewer's Shot): { ok, svg, widthMm, heightMm, marks, report, renders, warnings } */
    panel: (request) => ask('panel', { request: request?.scene ? { ...request, scene: absolute(request.scene) } : request }),
    /** the Stage's canvas size, CSS px: { ok, size: [w, h] } */
    stage: () => ask('stage', {}),
  }
}
