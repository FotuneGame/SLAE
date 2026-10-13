import { eliminate } from './index.js';

/** Определитель через исключение Гаусса (0 для вырожденной матрицы). */
export function det(A) {
  try {
    const { M, sign } = eliminate(A, false, 0);
    return M.reduce((p, row, i) => p * row[i], sign);
  } catch {
    return 0;
  }
}
