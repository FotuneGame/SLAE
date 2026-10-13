import { normMat, spectralRadius, iterationMatrix, seidelMatrix } from './index.js';

/** Матрицы итераций Якоби (B) и Зейделя (Bz), их нормы ‖·‖∞ и спектральные радиусы ρ. */
export function convergence(A) {
  const B = iterationMatrix(A);
  const Bz = seidelMatrix(B);
  return { B, Bz, rhoB: spectralRadius(B), rhoBz: spectralRadius(Bz), normB: normMat(B), normBz: normMat(Bz) };
}
