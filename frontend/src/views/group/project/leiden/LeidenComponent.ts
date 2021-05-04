import { ComponentContainer } from "golden-layout";
import Vue from "vue";
import LeidenView from "./LeidenView.vue";

export class LeidenComponent {
  static readonly typeName = LeidenComponent.name;

  private readonly _element: any;

  private readonly _containerClickListener = () => this.handleClickFocusEvent();
  private readonly _containerFocusinListener = () => this.handleClickFocusEvent();
  private readonly _beforeComponentReleaseEventListener = () => this.handleBeforeComponentReleaseEvent();
  private readonly _showEventListener = () => this.handleShowEvent();

  constructor(private _container: ComponentContainer, store: any, parent) {
    const ComponentClass = Vue.extend(LeidenView);
    this._element = new ComponentClass({ store: store, parent: parent });
    this._element.$mount();
    this._container.element.appendChild(this._element.$el);

    this._container.stateRequestEvent = () => this.handleContainerStateRequestEvent();
    this._container.addEventListener("beforeComponentRelease", this._beforeComponentReleaseEventListener);
    this._container.addEventListener("show", this._showEventListener);

    this._container.element.addEventListener("click", this._containerClickListener);
    this._container.element.addEventListener("focusin", this._containerFocusinListener);
  }

  private handleContainerStateRequestEvent(): string | undefined {
    return undefined;
  }

  private handleBeforeComponentReleaseEvent(): void {
    console.log("handleBeforeComponentReleaseEvent");
    this._element.$destroy();
    this._container.element.removeChild(this._element.$el);
    this._container.removeEventListener("show", this._showEventListener);
    this._container.removeEventListener("beforeComponentRelease", this._beforeComponentReleaseEventListener);
    this._container.element.removeEventListener("click", this._containerClickListener);
    this._container.element.removeEventListener("focusin", this._containerFocusinListener);
  }

  private handleShowEvent(): void {
    console.log("handleShowEvent");
  }

  private handleClickFocusEvent(): void {
    this._container.focus();
  }
}
