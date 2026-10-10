// График невязки в формате SVG (масштабы и линии строит d3).
import { range } from 'd3-array';
import { scaleLinear, scaleLog } from 'd3-scale';
import { line } from 'd3-shape';

const W = 720, H = 420, M = { l: 72, r: 24, t: 48, b: 84 };
const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const pow10 = (e) => '10' + [...String(e)].map((c) => SUP[c]).join('');

// series: [{ name, color, dash?, values: [невязка при k = 0, 1, 2, ...] }]
export function residualChart(title, series) {
  const all = series.flatMap((s) => s.values).map((v) => Math.max(v, 1e-16));
  const lo = Math.floor(Math.log10(Math.min(...all)));
  const hi = Math.max(lo + 1, Math.ceil(Math.log10(Math.max(...all))));
  const kMax = Math.max(1, ...series.map((s) => s.values.length - 1));
  const x = scaleLinear().domain([0, kMax]).range([M.l, W - M.r]);
  const y = scaleLog().domain([10 ** lo, 10 ** hi]).range([H - M.b, M.t]);
  const path = line().x((_, k) => x(k)).y((v) => y(Math.max(v, 1e-16)));
  const dash = (s) => (s.dash ? ` stroke-dasharray="${s.dash}"` : '');

  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="Arial, sans-serif" font-size="12">
        <rect width="${W}" height="${H}" fill="#fff"/>
        <text x="${W / 2}" y="26" text-anchor="middle" font-size="15" font-weight="bold">${title}</text>
        ${range(lo, hi + 1, Math.ceil((hi - lo) / 8)).map((e) => `<line x1="${M.l}" x2="${W - M.r}" y1="${y(10 ** e)}" y2="${y(10 ** e)}" stroke="#e5e7eb"/><text x="${M.l - 8}" y="${y(10 ** e) + 4}" text-anchor="end">${pow10(e)}</text>`).join('\n')}
        ${x.ticks(6).filter(Number.isInteger).map((k) => `<text x="${x(k)}" y="${H - M.b + 18}" text-anchor="middle">${k}</text>`).join('\n')}
        <rect x="${M.l}" y="${M.t}" width="${W - M.l - M.r}" height="${H - M.t - M.b}" fill="none" stroke="#9ca3af"/>
        <text x="${W / 2}" y="${H - 38}" text-anchor="middle">Номер итерации k</text>
        <text x="16" y="${H / 2}" text-anchor="middle" transform="rotate(-90 16 ${H / 2})">‖Ax − b‖∞</text>
        ${series.map((s, i) => `<path d="${path(s.values)}" fill="none" stroke="${s.color}" stroke-width="2"${dash(s)}/>
        <line x1="${W / 2 - 150 + i * 170}" x2="${W / 2 - 122 + i * 170}" y1="${H - 16}" y2="${H - 16}" stroke="${s.color}" stroke-width="2"${dash(s)}/>
        <text x="${W / 2 - 114 + i * 170}" y="${H - 12}">${s.name}</text>`).join('\n')}
        </svg>`;
}
