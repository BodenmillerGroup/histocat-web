import { Mutations } from "vuex-smart-module";
import { AnnotationsState } from ".";

export class AnnotationsMutations extends Mutations<AnnotationsState> {
  addCellClass(payload: { name: string; color: string }) {
    const newValue = { ...this.state.classes };
    newValue[payload.name] = payload.color;
    this.state.classes = newValue;
  }

  updateCellClass(payload: { prevName: string; nextName: string; color: string }) {
    const newValue = { ...this.state.classes };
    delete newValue[payload.prevName];
    newValue[payload.nextName] = payload.color;
    this.state.classes = newValue;
  }

  deleteCellClass(name: string) {
    const newValue = { ...this.state.classes };
    delete newValue[name];
    this.state.classes = newValue;
  }

  addAnnotation(payload: { name: string; cellClass: string; cells: string[] }) {
    this.state.annotations = this.state.annotations.concat({
      name: payload.name,
      cellClass: payload.cellClass,
      cells: payload.cells,
      visible: true,
    });
  }

  updateAnnotation(payload: { prevName: string; nextName: string; cellClass: string }) {
    const annotation = this.state.annotations.find((v) => v.name === payload.prevName);
    if (annotation) {
      annotation.name = payload.nextName;
      annotation.cellClass = payload.cellClass;
      this.state.annotations = this.state.annotations.slice();
    }
  }

  deleteAnnotation(name: string) {
    this.state.annotations = this.state.annotations.filter((v) => v.name !== name);
  }

  reset() {
    // acquire initial state
    const s = new AnnotationsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
