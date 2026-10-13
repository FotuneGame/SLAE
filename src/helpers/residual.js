import { dot, normVec } from './vector.js';

// Невязка ‖Ax − b‖∞
export const residual = (A, x, b) => normVec(A.map((row, i) => dot(row, x) - b[i]));
