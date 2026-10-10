// Практическая работа: решение СЛАУ, вариант 18. Запуск: node src/index.js
// Входные матрицы A и b задаются в файле input/data.json.
import { mkdirSync, writeFileSync } from 'node:fs';
import Table from 'cli-table3';
import { residualChart } from './chart.js';
import { systems, params, exactSolutions, exactSolutionProvided } from './data.js';
import {
  gauss, cramer, gaussJordan, inverse, normMat, normVec, residual,
  iterate, iterationMatrices, normalEquations,
} from './solvers.js';

const EPS = params.eps;
const STATUS = { converged: 'шаг < ε', diverged: 'расходится', maxIter: 'нет сходимости' };
const fmt = (v) => (v !== 0 && (Math.abs(v) < 1e-4 || Math.abs(v) >= 1e6) ? v.toExponential(2) : v.toFixed(4));
const table = (head, rows) => {
  const t = new Table({ head, colAligns: head.map(() => 'right'), style: { head: ['bold'], border: [] },
    chars: { mid: '', 'left-mid': '', 'mid-mid': '', 'right-mid': '' } });
  t.push(...rows.map((r) => r.map(String)));
  console.log(t.toString());
};

// Запуск итераций для системы (As, bs); невязка и ошибка считаются для исходной системы (A, b).
function run(A, b, exact, As, bs, opts) {
  const { xs, status } = iterate(As, bs, opts);
  const x = xs.at(-1);
  return { status, k: xs.length - 1, x, error: normVec(x.map((v, i) => v - exact[i])),
    residuals: xs.slice(0, 101).map((xk) => residual(A, xk, b)) };
}

const results = systems.map(({ A, b }, idx) => {
  console.log(`\n=== Система ${idx + 1} ===`);

  // Прямые методы
  const g = gauss(A, b), c = cramer(A, b), gj = gaussJordan(A, b);
  // Опорное решение x*: из data.json, а при его отсутствии — лучшее вычисленное (Гаусс–Жордан).
  const ref = exactSolutions[idx];
  if (!exactSolutionProvided) console.log(`x* не задано в data.json — ошибки считаются относительно лучшего вычисленного значения [${ref.map(fmt).join(', ')}]`);
  table(['метод', 'x1', 'x2', 'x3', 'x4', '‖Ax−b‖'], [['Гаусс', g.x], ['Крамер', c.x], ['Гаусс-Жордан', gj]]
    .map(([name, x]) => [name, ...x.map(fmt), fmt(residual(A, x, b))]));
  const Ainv = inverse(A);
  const cond = normMat(A) * normMat(Ainv);
  console.log(`det(A) = ${fmt(c.delta)}, ‖A‖ = ${fmt(normMat(A))}, ‖A⁻¹‖ = ${fmt(normMat(Ainv))}, cond(A) = ${cond.toFixed(2)}`);

  // Итерационные методы
  const m = iterationMatrices(A);
  const jac = run(A, b, ref, A, b, { eps: EPS });
  const sei = run(A, b, ref, A, b, { eps: EPS, seidel: true });
  const NE = normalEquations(A, b);
  const mN = iterationMatrices(NE.A);
  const bigRun = { eps: EPS, maxIter: 200000, blowUp: 1e30 };
  const jacN = run(A, b, ref, NE.A, NE.b, bigRun);
  const seiN = run(A, b, ref, NE.A, NE.b, { ...bigRun, seidel: true });

  console.log('\nИтерации (x⁰ = 0):');
  table(['метод', '‖B‖', 'ρ(B)', 'итераций', 'результат', '‖x−x*‖'], [
    ['Якоби', fmt(m.normB), fmt(m.rhoB), jac.k, STATUS[jac.status], fmt(jac.error)],
    ['Зейдель', fmt(m.normBz), fmt(m.rhoBz), sei.k, STATUS[sei.status], fmt(sei.error)],
    ['Якоби, AᵀAx=Aᵀb', '', fmt(mN.rhoB), jacN.k, STATUS[jacN.status], fmt(jacN.error)],
    ['Зейдель, AᵀAx=Aᵀb', '', fmt(mN.rhoBz), seiN.k, STATUS[seiN.status], fmt(seiN.error)],
  ]);
  console.log('\nНевязка ‖Ax(k)−b‖:');
  table(['k', 'Якоби', 'Зейдель'], [0, 1, 2, 3, 4, 5, 10].filter((k) => k < sei.residuals.length && k < jac.residuals.length)
    .map((k) => [k, fmt(jac.residuals[k]), fmt(sei.residuals[k])]));

  // Графики
  const series = (j, s) => [{ name: 'Метод Якоби', color: '#2563eb', values: j.residuals },
    { name: 'Метод Зейделя', color: '#dc2626', dash: '6 4', values: s.residuals }];
  mkdirSync('output', { recursive: true });
  writeFileSync(`output/system${idx + 1}.svg`, residualChart(`Система ${idx + 1}: исходная система`, series(jac, sei)));
  writeFileSync(`output/system${idx + 1}-AtA.svg`, residualChart(`Система ${idx + 1}: AᵀAx = Aᵀb (первые 100 итераций)`, series(jacN, seiN)));

  return { A, b, direct: { gauss: g.x, cramer: c.x, gaussJordan: gj, triangular: g.triangular, delta: c.delta, dets: c.dets,
      residuals: [g.x, c.x, gj].map((x) => residual(A, x, b)) },
    Ainv, normA: normMat(A), normAinv: normMat(Ainv), cond, matrices: m, jacobi: jac, seidel: sei,
    ata: { jacobi: jacN, seidel: seiN, rhoB: mN.rhoB, rhoBz: mN.rhoBz } };
});

writeFileSync('output/results.json', JSON.stringify(results, null, 2));
console.log('\nРезультаты и графики сохранены в папке output');
