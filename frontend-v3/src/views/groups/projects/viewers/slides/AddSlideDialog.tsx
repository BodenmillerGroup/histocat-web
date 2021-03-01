import {
  Button,
  Classes,
  Dialog,
  FileInput,
  FormGroup,
  Intent,
} from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { useState } from "react";
import { useProjectsStore } from "modules/projects";

type AddSlideDialogProps = {
  isOpen: boolean;
  handleClose(): void;
};

export function AddSlideDialog(props: AddSlideDialogProps) {
  const { handleSubmit } = useForm();
  const uploadSlide = useProjectsStore((state) => state.uploadSlide);
  const [file, setFile] = useState<File | null>(null);

  const onSubmit = async (values: any) => {
    if (file) {
      // form is valid
      const formData = new FormData();
      formData.append("file", file);
      await uploadSlide(formData);
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
      title="Upload Slide"
      usePortal={true}
      isOpen={props.isOpen}
      className={Classes.DARK}
      canOutsideClickClose={false}
    >
      <form onSubmit={handleSubmit(onSubmit)}>
        <div className={Classes.DIALOG_BODY}>
          <FormGroup label="Slide file" labelFor="file-input">
            <FileInput
              id="file-input"
              text={file ? file.name : "Choose slide file..."}
              fill={true}
              onInputChange={onFileChange}
              hasSelection={!!file}
              inputProps={{
                accept: ".zip,.mcd",
              }}
            />
          </FormGroup>
        </div>
        <div className={Classes.DIALOG_FOOTER}>
          <div className={Classes.DIALOG_FOOTER_ACTIONS}>
            <Button onClick={props.handleClose} text="Cancel" />
            <Button type="reset" text="Reset" onClick={() => setFile(null)} />
            <Button type="submit" intent={Intent.PRIMARY} text="Upload" disabled={!file} />
          </div>
        </div>
      </form>
    </Dialog>
  );
}
