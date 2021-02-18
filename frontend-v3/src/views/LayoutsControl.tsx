import { Button, ControlGroup, HTMLSelect, Intent } from "@blueprintjs/core";
import shallow from "zustand/shallow";
import { useLayoutsStore } from "../modules/layouts";

export function LayoutsControl() {
  const { layouts, getLayout, addLayout, deleteLayout } = useLayoutsStore(
    (state) => ({
      layouts: state.layouts,
      getLayout: state.getLayout,
      addLayout: state.addLayout,
      deleteLayout: state.deleteLayout,
    }),
    shallow
  );

  return (
    <ControlGroup vertical={false}>
      <HTMLSelect>
        {layouts.map((layout) => {
          return <option>{layout.name}</option>;
        })}
      </HTMLSelect>
      <Button text="Load layout" icon="layout" intent={Intent.PRIMARY} />
      <Button icon="delete" intent={Intent.DANGER} />
    </ControlGroup>
  );
}
