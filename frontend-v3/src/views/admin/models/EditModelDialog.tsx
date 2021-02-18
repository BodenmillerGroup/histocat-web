import { Button, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { useModelsStore } from "modules/models";
import { IModel, IModelUpdate } from "modules/models/models";

type EditModelDialogProps = {
  isOpen: boolean;
  handleClose(): void;
  model: IModel;
};

export function EditModelDialog(props: EditModelDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const updateModel = useModelsStore((state) => state.updateModel);

  const onSubmit = async (values: any) => {
    // form is valid
    const params: IModelUpdate = {
      name: values.name,
      description: values.description,
    };
    await updateModel(props.model.id, params);
    props.handleClose();
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Edit Model"
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
              placeholder="Enter model name"
              defaultValue={props.model.name}
              inputRef={register({
                required: "Name is required",
              })}
            />
          </FormGroup>

          <FormGroup label="Description" labelFor="description-input">
            <InputGroup
              id="description-input"
              name="description"
              placeholder="Enter model description"
              defaultValue={props.model.description}
              inputRef={register({})}
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
