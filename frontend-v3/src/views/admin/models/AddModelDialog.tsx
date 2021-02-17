import { Button, Classes, Dialog, FileInput, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { useModelsStore } from "modules/models";
import { useState } from "react";

type AddModelDialogProps = {
  isOpen: boolean;
  handleClose(): void;
};

export function AddModelDialog(props: AddModelDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const createModel = useModelsStore((state) => state.createModel);
  const [file, setFile] = useState<File | null>(null);

  const onSubmit = async (values: any) => {
    if (file) {
      // form is valid
      const formData = new FormData();
      formData.append("name", values.name);
      formData.append("description", values.description);
      formData.append("file", file);
      await createModel(formData);
      setFile(null);
      props.handleClose();
    }
  };

  const onFileChange = (event: React.FormEvent<HTMLInputElement>) => {
    if ((event.target as HTMLInputElement).files !== null) {
      const file = (event.target as HTMLInputElement).files![0];
      setFile(file);
    }
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Add Model"
      usePortal={true}
      isOpen={props.isOpen}
      className="bp3-dark"
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
              inputRef={register({})}
            />
          </FormGroup>

          <FormGroup label="Model file" labelFor="file-input">
            <FileInput
              id="file-input"
              text={file ? file.name : "Choose model file..."}
              fill={true}
              onInputChange={onFileChange}
              hasSelection={!!file}
              inputProps={{
                accept: ".zip",
              }}
            />
          </FormGroup>
        </div>
        <div className={Classes.DIALOG_FOOTER}>
          <div className={Classes.DIALOG_FOOTER_ACTIONS}>
            <Button onClick={props.handleClose} text="Cancel" />
            <Button type="reset" text="Reset" onClick={() => setFile(null)} />
            <Button type="submit" intent={Intent.PRIMARY} text="Save" />
          </div>
        </div>
      </form>
    </Dialog>
  );
}
