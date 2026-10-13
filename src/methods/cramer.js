import { det } from '../helpers/index.js';

/** Метод Крамера: xᵢ = det(Aᵢ) / det(A). */
export function cramer(A, b) {
  const delta = det(A);
  if (Math.abs(delta) < 1e-14) throw new Error('det(A) = 0, метод Крамера неприменим');
  const dets = A.map((_, i) => det(A.map((row, k) => row.map((a, j) => (j === i ? b[k] : a)))));
  return { x: dets.map((d) => d / delta), delta, dets };
}
