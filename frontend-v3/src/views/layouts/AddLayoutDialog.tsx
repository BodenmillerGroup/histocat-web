import { Button, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { useLayoutsStore } from "modules/layouts";

type AddLayoutDialogProps = {
  isOpen: boolean;
  handleClose(): void;
};

export function AddLayoutDialog(props: AddLayoutDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const addLayout = useLayoutsStore(state => state.addLayout);

  const onSubmit = async (values: any) => {
    // form is valid
    addLayout(values.name);
    props.handleClose();
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Add Layout"
      usePortal={true}
      isOpen={props.isOpen}
      className={Classes.DARK}
      canOutsideClickClose={false}
    >
      <form onSubmit={handleSubmit(onSubmit)}>
        <div className={Classes.DIALOG_BODY}>
          <FormGroup
            label="Name"
            labelFor="name-input"
            labelInfo="(required)"
            intent="danger"
            helperText={errors?.name?.message}
          >
            <InputGroup
              id="name-input"
              name="name"
              placeholder="Enter name"
              inputRef={register({
                required: "Name is required",
              })}
            />
          </FormGroup>
        </div>
        <div className={Classes.DIALOG_FOOTER}>
          <div className={Classes.DIALOG_FOOTER_ACTIONS}>
            <Button onClick={props.handleClose} text="Cancel" />
            <Button type="reset" text="Reset" />
            <Button type="submit" intent={Intent.PRIMARY} text="Save" />
          </div>
        </div>
      </form>
    </Dialog>
  );
}
