// Исследование результатов работы алгоритмов при случайном изменении матрицы A и/или
// вектора b в пределах ±percent%. Прогонов много (runs), изменения случайные — метод Монте-Карло.
// Запуск: npm run start:random
import { mkdirSync, writeFileSync } from 'node:fs';
import Table from 'cli-table3';
import { systems, exactSolutions, exactSolutionProvided, params } from './data.js';
import {
  gauss, cramer, gaussJordan, inverse, normMat, normVec, residual,
  iterate, iterationMatrices, normalEquations,
} from './solvers.js';

const { percent, runs, seed, maxIter } = params.random;
const EPS = params.eps;

const STATUS = { converged: 'шаг < ε', diverged: 'расходится', maxIter: 'нет сходимости' };
const fmt = (v) => (!Number.isFinite(v) ? '—' : v !== 0 && (Math.abs(v) < 1e-4 || Math.abs(v) >= 1e6) ? v.toExponential(2) : v.toFixed(4));
const table = (head, rows) => {
  const t = new Table({ head, colAligns: head.map(() => 'right'), style: { head: ['bold'], border: [] },
    chars: { mid: '', 'left-mid': '', 'mid-mid': '', 'right-mid': '' } });
  t.push(...rows.map((r) => r.map(String)));
  console.log(t.toString());
};

// Детерминированный ГПСЧ (mulberry32): прогоны воспроизводимы при одном seed.
function mulberry32(a) {
  return () => {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Каждый элемент умножается на (1 + δ), δ ∈ [−p/100, p/100] — относительное изменение ±p%.
const perturbMat = (M, p, rnd) => M.map((row) => row.map((v) => v * (1 + (rnd() * 2 - 1) * p / 100)));
const perturbVec = (v, p, rnd) => v.map((x) => x * (1 + (rnd() * 2 - 1) * p / 100));

const stats = (arr) => {
  if (!arr.length) return { mean: NaN, min: NaN, max: NaN, std: NaN, n: 0 };
  const mean = arr.reduce((s, v) => s + v, 0) / arr.length;
  return { mean, min: Math.min(...arr), max: Math.max(...arr), n: arr.length,
    std: Math.sqrt(arr.reduce((s, v) => s + (v - mean) ** 2, 0) / arr.length) };
};

// «Изменять по очереди и вместе»: три случая варьирования.
const CASES = [
  { name: 'только A', changeA: true, changeB: false },
  { name: 'только b', changeA: false, changeB: true },
  { name: 'A и b вместе', changeA: true, changeB: true },
];

// Один прогон: решаем систему всеми методами для (возможно изменённых) A, b.
function analyze(A, b) {
  const g = gauss(A, b), c = cramer(A, b), gj = gaussJordan(A, b);
  const cond = normMat(A) * normMat(inverse(A));
  const m = iterationMatrices(A);
  const jac = iterate(A, b, { eps: EPS });
  const sei = iterate(A, b, { eps: EPS, seidel: true });
  const NE = normalEquations(A, b);
  const jacN = iterate(NE.A, NE.b, { eps: EPS, maxIter, blowUp: 1e30 });
  const seiN = iterate(NE.A, NE.b, { eps: EPS, seidel: true, maxIter, blowUp: 1e30 });
  return { g, c, gj, cond, m, jac, sei, jacN, seiN };
}

// Накопители статистики по одному случаю.
const newAcc = () => ({
  failed: 0,
  cond: [],
  rhoB: [], rhoBz: [],
  direct: { gauss: [], cramer: [], gaussJordan: [] },  // ошибки ‖x−x*‖
  res: { gauss: [], cramer: [], gaussJordan: [] },     // невязки ‖Ax−b‖
  iter: {
    jacobi: { err: [], iter: [], conv: 0, total: 0 },
    seidel: { err: [], iter: [], conv: 0, total: 0 },
    jacobiAtA: { err: [], iter: [], conv: 0, total: 0 },
    seidelAtA: { err: [], iter: [], conv: 0, total: 0 },
  },
});

// Ошибка ‖x − x*‖ относительно опорного решения ref (заданного или вычисленного).
const errOf = (x, ref) => normVec(x.map((v, i) => v - ref[i]));
const recordIter = (bucket, xs, status, ref) => {
  bucket.total++;
  if (status === 'converged') {
    bucket.conv++;
    bucket.iter.push(xs.length - 1);
    bucket.err.push(errOf(xs.at(-1), ref));
  }
};

const rnd = mulberry32(seed);
const results = systems.map(({ A: A0, b: b0 }, idx) => {
  // Опорное решение: из data.json либо лучшее вычисленное (Гаусс–Жордан исходной системы).
  const ref = exactSolutions[idx];

  console.log(`\n=== Система ${idx + 1}: метод Монте-Карло, ${runs} прогонов, изменение ±${percent}% (seed = ${seed}) ===`);
  if (!exactSolutionProvided) console.log(`x* не задано в data.json — используется лучшее вычисленное значение: [${ref.map(fmt).join(', ')}]`);
  console.log(`Итерационные методы: ε = ${EPS}, не более ${maxIter} итераций на прогон.`);

  const cases = CASES.map((cs) => {
    const acc = newAcc();
    for (let r = 0; r < runs; r++) {
      const A = cs.changeA ? perturbMat(A0, percent, rnd) : A0;
      const b = cs.changeB ? perturbVec(b0, percent, rnd) : b0;
      let res;
      try { res = analyze(A, b); } catch { acc.failed++; continue; }

      acc.cond.push(res.cond);
      acc.rhoB.push(res.m.rhoB);
      acc.rhoBz.push(res.m.rhoBz);
      for (const [key, x] of [['gauss', res.g.x], ['cramer', res.c.x], ['gaussJordan', res.gj]]) {
        acc.direct[key].push(errOf(x, ref));
        acc.res[key].push(residual(A, x, b));
      }
      recordIter(acc.iter.jacobi, res.jac.xs, res.jac.status, ref);
      recordIter(acc.iter.seidel, res.sei.xs, res.sei.status, ref);
      recordIter(acc.iter.jacobiAtA, res.jacN.xs, res.jacN.status, ref);
      recordIter(acc.iter.seidelAtA, res.seiN.xs, res.seiN.status, ref);
    }
    return { ...cs, acc };
  });

  // Таблица 1: обусловленность и прямые методы. В каждой ячейке — мин. / сред. / макс. за прогон.
  console.log('\nПрямые методы: ‖x − x*‖ и невязка ‖Ax − b‖ (мин. / сред. / макс. за прогон):');
  const triple = (s) => `${fmt(s.min)} / ${fmt(s.mean)} / ${fmt(s.max)}`;
  table(['изменение', 'cond(A), мин/сред/макс', 'Гаусс ‖x−x*‖, мин/сред/макс', 'Крамер ‖x−x*‖, мин/сред/макс',
    'Г-Жордан ‖x−x*‖, мин/сред/макс', 'невязка, мин/сред/макс'],
    cases.map(({ name, acc }) => {
      const n = acc.res.gauss.length;
      const resWorst = Array.from({ length: n }, (_, i) => Math.max(acc.res.gauss[i], acc.res.cramer[i], acc.res.gaussJordan[i]));
      return [name, triple(stats(acc.cond)), triple(stats(acc.direct.gauss)), triple(stats(acc.direct.cramer)),
        triple(stats(acc.direct.gaussJordan)), triple(stats(resWorst))];
    }));

  // Таблица 2: итерационные методы. В каждой ячейке — мин. / сред. / макс. по сошедшимся прогонам.
  console.log('\nИтерационные методы (сходимость, итерации и ‖x − x*‖ по сошедшимся прогонам):');
  const iterRows = (acc, name) => {
    const labels = [['jacobi', 'Якоби'], ['seidel', 'Зейдель'], ['jacobiAtA', 'Якоби, AᵀAx=Aᵀb'], ['seidelAtA', 'Зейдель, AᵀAx=Aᵀb']];
    return labels.map(([key, label]) => {
      const bk = acc.iter[key];
      const pct = ((bk.conv / bk.total) * 100).toFixed(0);
      return [name, label, `${pct}%`, `${bk.conv}/${bk.total}`, triple(stats(bk.iter)), triple(stats(bk.err))];
    });
  };
  table(['изменение', 'метод', 'сходимость', 'сошлось', 'итер., мин/сред/макс', '‖x−x*‖, мин/сред/макс'],
    cases.flatMap(({ name, acc }) => iterRows(acc, name)));
  console.log(`ρ(B) Якоби: ${cases.map((c) => fmt(stats(c.acc.rhoB).mean)).join(' / ')} (только A / только b / A и b); ` +
    `ρ(Bz) Зейделя: ${cases.map((c) => fmt(stats(c.acc.rhoBz).mean)).join(' / ')}`);

  return { percent, runs, seed, cases: cases.map(({ name, changeA, changeB, acc }) => ({
    case: name, changeA, changeB, failed: acc.failed,
    cond: stats(acc.cond), rhoB: stats(acc.rhoB), rhoBz: stats(acc.rhoBz),
    direct: Object.fromEntries(Object.entries(acc.direct).map(([k, v]) => [k, stats(v)])),
    residual: Object.fromEntries(Object.entries(acc.res).map(([k, v]) => [k, stats(v)])),
    iter: Object.fromEntries(Object.entries(acc.iter).map(([k, v]) => [k, {
      conv: v.conv, total: v.total, iter: stats(v.iter), err: stats(v.err) }])),
  })) };
});

mkdirSync('output', { recursive: true });
writeFileSync('output/random-results.json', JSON.stringify(results, null, 2));
console.log('\nСтатистика случайных прогонов сохранена в output/random-results.json');
