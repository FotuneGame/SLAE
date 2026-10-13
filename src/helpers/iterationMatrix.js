import { matMul, inverse } from './index.js';

/** Матрица B = −D⁻¹(L+U) итерационного метода Якоби: bᵢⱼ = −aᵢⱼ/aᵢᵢ (i≠j), bᵢᵢ = 0. */
export const iterationMatrix = (A) => A.map((row, i) => row.map((a, j) => (i === j ? 0 : -a / row[i])));

/** Матрица Bz = (I − L_B)⁻¹ U_B метода Зейделя, где B = L_B + U_B. */
export function seidelMatrix(B) {
  const part = (keep) => B.map((row, i) => row.map((v, j) => (keep(i, j) ? v : 0)));
  const IminusL = part((i, j) => j < i).map((row, i) => row.map((v, j) => +(i === j) - v));
  return matMul(inverse(IminusL), part((i, j) => j > i));
}
