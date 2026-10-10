// Решение СЛАУ Ax = b. Матрица — массив строк, вектор — обычный массив.

const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
const transpose = (A) => A[0].map((_, j) => A.map((row) => row[j]));
const matMul = (A, B) => A.map((row) => transpose(B).map((col) => dot(row, col)));
const augment = (A, b) => A.map((row, i) => [...row, ...(Array.isArray(b[i]) ? b[i] : [b[i]])]);

export const normVec = (v) => Math.max(...v.map(Math.abs));                                   // ‖v‖∞
export const normMat = (A) => Math.max(...A.map((row) => row.reduce((s, a) => s + Math.abs(a), 0))); // ‖A‖∞
export const residual = (A, x, b) => normVec(A.map((row, i) => dot(row, x) - b[i]));           // ‖Ax−b‖∞

// ================= Прямые методы =================

// Прямой ход с выбором главного элемента по столбцу. jordan=true — ещё и обратный ход (Гаусс–Жордан).
function eliminate(M, jordan = false, tol = 1e-14) {
  const n = M.length;
  let sign = 1;
  M = M.map((row) => [...row]);
  for (let k = 0; k < n; k++) {
    const p = M.map((row, i) => (i < k ? -1 : Math.abs(row[k]))).reduce((best, v, i, a) => (v > a[best] ? i : best), k);
    if (Math.abs(M[p][k]) <= tol) throw new Error('Матрица вырождена');
    if (p !== k) { [M[p], M[k]] = [M[k], M[p]]; sign = -sign; }
    if (jordan) M[k] = M[k].map((v) => v / M[k][k]);              // делим строку на ведущий элемент
    M.forEach((row, i) => {
      if (i === k || (!jordan && i < k)) return;
      const f = row[k] / M[k][k];                                   // обнуляем столбец k
      M[i] = row.map((v, j) => v - f * M[k][j]);
    });
  }
  return { M, sign };
}

export function gauss(A, b) {
  const n = A.length;
  const { M } = eliminate(augment(A, b));                           // верхнетреугольный вид
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {                                // обратный ход
    x[i] = (M[i][n] - dot(M[i].slice(i + 1, n), x.slice(i + 1))) / M[i][i];
  }
  return { x, triangular: M };
}

export const gaussJordan = (A, b) => eliminate(augment(A, b), true).M.map((row) => row[row.length - 1]);

export function det(A) {
  try {
    const { M, sign } = eliminate(A, false, 0);
    return M.reduce((p, row, i) => p * row[i], sign);
  } catch {
    return 0;
  }
}

export function cramer(A, b) {
  const delta = det(A);
  if (Math.abs(delta) < 1e-14) throw new Error('det(A) = 0, метод Крамера неприменим');
  const dets = A.map((_, i) => det(A.map((row, k) => row.map((a, j) => (j === i ? b[k] : a)))));
  return { x: dets.map((d) => d / delta), delta, dets };
}

// Обратная матрица: [A | I] → [I | A⁻¹]
export const inverse = (A) => {
  const I = A.map((_, i) => A.map((_, j) => +(i === j)));
  return eliminate(augment(A, I), true).M.map((row) => row.slice(A.length));
};

// ================= Итерационные методы =================

// Итерации по Якоби (seidel = false) или Зейделю (seidel = true) из x⁽⁰⁾ = 0.
// Возвращает все приближения xs[k] и причину остановки.
export function iterate(A, b, { seidel = false, eps = 1e-6, maxIter = 1000, blowUp = 1e12 } = {}) {
  let x = new Array(A.length).fill(0);
  const xs = [x];
  for (let k = 1; k <= maxIter; k++) {
    const next = [...x];
    A.forEach((row, i) => {
      const old = seidel ? next : x;                                // Зейдель берёт уже обновлённые значения
      next[i] = (b[i] - dot(row, old) + row[i] * old[i]) / row[i];  // xᵢ = (bᵢ − Σ_{j≠i} aᵢⱼxⱼ) / aᵢᵢ
    });
    xs.push(next);
    const step = normVec(next.map((v, i) => v - x[i]));
    x = next;
    if (!(normVec(x) < blowUp)) return { xs, status: 'diverged' };
    if (step < eps) return { xs, status: 'converged' };
  }
  return { xs, status: 'maxIter' };
}

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

// Матрица B = −D⁻¹(L+U) метода Якоби и матрица Bz = (I − L_B)⁻¹ U_B метода Зейделя (B = L_B + U_B).
export function iterationMatrices(A) {
  const B = A.map((row, i) => row.map((a, j) => (i === j ? 0 : -a / row[i])));
  const part = (keep) => B.map((row, i) => row.map((v, j) => (keep(i, j) ? v : 0)));
  const IminusL = part((i, j) => j < i).map((row, i) => row.map((v, j) => +(i === j) - v));
  const Bz = matMul(inverse(IminusL), part((i, j) => j > i));
  return { B, Bz, rhoB: spectralRadius(B), rhoBz: spectralRadius(Bz), normB: normMat(B), normBz: normMat(Bz) };
}

// Система AᵀAx = Aᵀb: симметричная и положительно определённая, для неё Зейдель сходится.
export const normalEquations = (A, b) => ({ A: matMul(transpose(A), A), b: transpose(A).map((col) => dot(col, b)) });
