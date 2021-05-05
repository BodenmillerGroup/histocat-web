import { ComponentContainer } from "golden-layout";

export abstract class LayoutComponent {
  protected readonly _component: any;

  private readonly _containerClickListener = () => this.handleClickFocusEvent();
  private readonly _containerFocusinListener = () => this.handleClickFocusEvent();
  private readonly _beforeComponentReleaseEventListener = () => this.handleBeforeComponentReleaseEvent();
  private readonly _showEventListener = () => this.handleShowEvent();
  private readonly _resizeEventListener = () => this.handleResizeEvent();

  constructor(protected _container: ComponentContainer, ComponentClass: any, store: any, parent) {
    this._component = new ComponentClass({ store: store, parent: parent });
    this._component.$mount();
    this._container.element.appendChild(this._component.$el);

    this._container.stateRequestEvent = () => this.handleContainerStateRequestEvent();
    this._container.addEventListener("beforeComponentRelease", this._beforeComponentReleaseEventListener);
    this._container.addEventListener("show", this._showEventListener);
    this._container.addEventListener("resize", this._resizeEventListener);

    this._container.element.addEventListener("click", this._containerClickListener);
    this._container.element.addEventListener("focusin", this._containerFocusinListener);
  }

  protected handleContainerStateRequestEvent(): string | undefined {
    return undefined;
  }

  protected handleBeforeComponentReleaseEvent(): void {
    console.log("handleBeforeComponentReleaseEvent");
    this._component.$destroy();
    this._container.element.removeChild(this._component.$el);
    this._container.removeEventListener("show", this._showEventListener);
    this._container.removeEventListener("resize", this._resizeEventListener);
    this._container.removeEventListener("beforeComponentRelease", this._beforeComponentReleaseEventListener);
    this._container.element.removeEventListener("click", this._containerClickListener);
    this._container.element.removeEventListener("focusin", this._containerFocusinListener);
  }

  protected handleShowEvent(): void {
    console.log("handleShowEvent");
  }

  protected handleResizeEvent(): void {
    console.log("handleResizeEvent");
  }

  protected handleClickFocusEvent(): void {
    this._container.focus();
  }
}
