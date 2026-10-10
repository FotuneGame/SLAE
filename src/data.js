// Загрузка входных данных варианта из input/data.json.
import { readFileSync } from 'node:fs';
import { gaussJordan } from './solvers.js';

const raw = JSON.parse(readFileSync(new URL('../input/data.json', import.meta.url), 'utf8'));

export const params = raw.params;
export const systems = raw.systems.map(({ A, b }) => ({ A, b }));

// Точное решение x*. Если оно не задано в data.json, вместо него берётся лучшее
// вычисленное значение — решение исходной системы прямым методом Гаусса–Жордана.
export const exactSolutionProvided = Array.isArray(raw.exactSolution);
export const exactSolutions = systems.map(({ A, b }) =>
  exactSolutionProvided ? raw.exactSolution : gaussJordan(A, b));
