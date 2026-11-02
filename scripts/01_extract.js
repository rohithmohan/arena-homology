#!/usr/bin/env node
// Extract the three embedded data literals (Cd = targets, He = scores, Pd = models)
// from the CoFold Arena JS bundle. Usage:
//   node 01_extract.js <bundle.js> <outdir>
// The bundle is untrusted: we only brace-match the literal text and evaluate that
// substring inside an empty vm context (no require / process / globals).
const fs = require('fs');
const vm = require('vm');
const path = require('path');

const [bundlePath, outDir] = process.argv.slice(2);
const src = fs.readFileSync(bundlePath, 'utf8');

function extractLiteral(marker) {
  const idx = src.indexOf(marker);
  if (idx < 0) throw new Error(`marker not found: ${marker}`);
  const start = idx + marker.indexOf('=') + 1;
  const open = src[start];
  const close = open === '[' ? ']' : '}';
  if (open !== '[' && open !== '{') throw new Error(`unexpected opener ${open} for ${marker}`);
  let depth = 0, inStr = null, esc = false;
  for (let i = start; i < src.length; i++) {
    const c = src[i];
    if (inStr) {
      if (esc) esc = false;
      else if (c === '\\') esc = true;
      else if (c === inStr) inStr = null;
      continue;
    }
    if (c === '"' || c === "'" || c === '`') { inStr = c; continue; }
    if (c === '[' || c === '{') depth++;
    else if (c === ']' || c === '}') {
      depth--;
      if (depth === 0) {
        if (c !== close) throw new Error('bracket mismatch');
        return src.slice(start, i + 1);
      }
    }
  }
  throw new Error(`unterminated literal for ${marker}`);
}

function evalLiteral(text) {
  // Literal must only contain data: reject identifiers that look like calls/functions.
  if (/=>|function\s*\(|\bnew\s+[A-Z]/.test(text)) throw new Error('literal contains code');
  return vm.runInNewContext('(' + text + ')', Object.create(null), { timeout: 5000 });
}

const targets = evalLiteral(extractLiteral('Cd=[{pdb'));
const scores = evalLiteral(extractLiteral('He={"'));
const models = evalLiteral(extractLiteral('Pd=[{id'));

fs.mkdirSync(outDir, { recursive: true });
fs.writeFileSync(path.join(outDir, 'targets.json'), JSON.stringify(targets, null, 1));
fs.writeFileSync(path.join(outDir, 'scores.json'), JSON.stringify(scores));
fs.writeFileSync(path.join(outDir, 'models.json'), JSON.stringify(models, null, 1));
console.log(`targets=${targets.length} scoreKeys=${Object.keys(scores).length} models=${models.length}`);
