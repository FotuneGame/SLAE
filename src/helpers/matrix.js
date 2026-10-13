// Операции над матрицами (матрица — массив строк).
import { dot } from './vector.js';

export const transpose = (A) => A[0].map((_, j) => A.map((row) => row[j]));
export const matMul = (A, B) => A.map((row) => transpose(B).map((col) => dot(row, col)));
export const augment = (A, b) => A.map((row, i) => [...row, ...(Array.isArray(b[i]) ? b[i] : [b[i]])]); // [A | b]
export const identity = (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => +(i === j)));
export const normMat = (A) => Math.max(...A.map((row) => row.reduce((s, a) => s + Math.abs(a), 0))); // ‖A‖∞
