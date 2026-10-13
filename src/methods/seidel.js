import { dot, normVec } from '../helpers/index.js';

/**
 * Метод Зейделя из x⁽⁰⁾ = 0: при вычислении xᵢ используются уже обновлённые xⱼ⁽ᵏ⁺¹⁾, j < i.
 * Возвращает все приближения xs[k] и причину остановки: converged | diverged | maxIter.
 */
export function seidel(A, b, { eps = 1e-6, maxIter = 1000, blowUp = 1e12 } = {}) {
  let x = new Array(A.length).fill(0);
  const xs = [x];
  for (let k = 1; k <= maxIter; k++) {
    const next = [...x];
    A.forEach((row, i) => {
      next[i] = (b[i] - dot(row, next) + row[i] * next[i]) / row[i];
    });
    xs.push(next);
    const step = normVec(next.map((v, i) => v - x[i]));
    x = next;
    if (!(normVec(x) < blowUp)) return { xs, status: 'diverged' };
    if (step < eps) return { xs, status: 'converged' };
  }
  return { xs, status: 'maxIter' };
}
