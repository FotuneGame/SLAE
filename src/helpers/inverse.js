import { augment, eliminate, identity } from './index.js';

/** Обратная матрица: [A | I] → [I | A⁻¹]. */
export const inverse = (A) => eliminate(augment(A, identity(A.length)), true).M.map((row) => row.slice(A.length));
