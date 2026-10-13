import { augment, eliminate } from '../helpers/index.js';

/** Метод Гаусса–Жордана: [A | b] → [E | x]. */
export const gaussJordan = (A, b) => eliminate(augment(A, b), true).M.map((row) => row[row.length - 1]);
