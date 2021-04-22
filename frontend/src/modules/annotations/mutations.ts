import { Mutations } from "vuex-smart-module";
import { AnnotationsState, defaultCellClasses } from ".";
import { IAnnotation } from "@/modules/annotations/models";

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

  setCellClasses(cellClasses: { [name: string]: string }) {
    this.state.cellClasses = cellClasses;
  }

  resetCellClasses() {
    this.state.cellClasses = { ...defaultCellClasses };
  }

  setAnnotations(annotations: IAnnotation[]) {
    this.state.annotations = annotations;
  }

  addAnnotation(payload: { cellClass: string; cellIds: string[] }) {
    const cellSet = new Set(payload.cellIds);
    const newAnnotations = [...this.state.annotations];
    // Calculate sets difference
    newAnnotations.forEach((annotation) => {
      annotation.cellIds = annotation.cellIds.filter((x) => !cellSet.has(x));
    });
    this.state.annotations = newAnnotations.concat({
      cellClass: payload.cellClass,
      cellIds: Array.from(cellSet),
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
