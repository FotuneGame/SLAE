// Операции над векторами (вектор — обычный массив).

export const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);         // (a, b)
export const normVec = (v) => Math.max(...v.map(Math.abs));                  // ‖v‖∞
