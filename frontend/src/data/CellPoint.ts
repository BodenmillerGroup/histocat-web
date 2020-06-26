export class CellPoint {
  constructor(
    public readonly acquisitionId: number,
    public readonly cellId: number,
    public readonly objectNumber: number,
    public readonly x: number,
    public readonly y: number,
    public readonly value: number
  ) {}
}
