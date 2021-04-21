import { Mutations } from "vuex-smart-module";
import { AnnotationsState, defaultCellClasses } from ".";

export class AnnotationsMutations extends Mutations<AnnotationsState> {
  addCellClass(payload: { name: string; color: string }) {
    const newValue = { ...this.state.cellClasses };
    newValue[payload.name] = payload.color;
    this.state.cellClasses = newValue;
  }

  updateCellClass(payload: { prevName: string; nextName: string; color: string }) {
    const newValue = { ...this.state.cellClasses };
    delete newValue[payload.prevName];
    newValue[payload.nextName] = payload.color;
    this.state.cellClasses = newValue;
  }

  deleteCellClass(name: string) {
    const newValue = { ...this.state.cellClasses };
    delete newValue[name];
    this.state.cellClasses = newValue;
  }

  resetCellClasses() {
    this.state.cellClasses = { ...defaultCellClasses };
  }

  addAnnotation(payload: { cellClass: string; cells: string[] }) {
    const cellSet = new Set(payload.cells);
    const newAnnotations = [...this.state.annotations];
    // Calculate sets difference
    newAnnotations.forEach((annotation) => {
      annotation.cells = new Set([...annotation.cells].filter((x) => !cellSet.has(x)));
    });
    this.state.annotations = newAnnotations.concat({
      cellClass: payload.cellClass,
      cells: cellSet,
      visible: true,
    });
  }

  updateAnnotation(payload: { index: number; cellClass: string }) {
    const annotation = this.state.annotations[payload.index];
    if (annotation) {
      annotation.cellClass = payload.cellClass;
      this.state.annotations = this.state.annotations.slice();
    }
  }

  deleteAnnotation(index: number) {
    const newValue = [...this.state.annotations];
    newValue.splice(index, 1);
    this.state.annotations = newValue;
  }

  reset() {
    // acquire initial state
    const s = new AnnotationsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
