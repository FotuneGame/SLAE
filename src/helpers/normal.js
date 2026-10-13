import { dot, matMul, transpose } from './index.js';

/** Нормальная система AᵀAx = Aᵀb: симметричная и положительно определённая, для неё Зейдель сходится. */
export const normal = (A, b) => ({ A: matMul(transpose(A), A), b: transpose(A).map((col) => dot(col, b)) });
