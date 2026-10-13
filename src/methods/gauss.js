import { augment, dot, eliminate } from '../helpers/index.js';

/** Метод Гаусса: прямой ход с выбором главного элемента + обратный ход. */
export function gauss(A, b) {
  const n = A.length;
  const { M } = eliminate(augment(A, b));                           // верхнетреугольный вид
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {                                // обратный ход
    x[i] = (M[i][n] - dot(M[i].slice(i + 1, n), x.slice(i + 1))) / M[i][i];
  }
  return { x, triangular: M };
}
