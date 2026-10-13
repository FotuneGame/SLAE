// Прямой ход с выбором главного элемента по столбцу. jordan=true — ещё и обратный ход (Гаусс–Жордан).
// Общая процедура для методов Гаусса, Гаусса–Жордана, определителя и обратной матрицы.
export function eliminate(M, jordan = false, tol = 1e-14) {
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
