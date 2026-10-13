import { matMul, normMat } from './matrix.js';

// Спектральный радиус: ρ(B) = lim ‖B^k‖^(1/k). Возводим B в степень 2^40 повторным возведением в квадрат
// с нормировкой (логарифм масштаба копится в logScale).
export function spectralRadius(B) {
  let M = B, logScale = 0;
  for (let j = 1; j <= 40; j++) {
    M = matMul(M, M);
    const s = normMat(M);
    M = M.map((row) => row.map((v) => v / s));
    logScale = 2 * logScale + Math.log(s);
  }
  return Math.exp(logScale / 2 ** 40);
}
