import { Button, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { useGroupsStore } from "modules/groups";
import { IGroup, IGroupUpdate } from "modules/groups/models";

type EditModelDialogProps = {
  isOpen: boolean;
  handleClose(): void;
  group: IGroup;
};

export function EditGroupDialog(props: EditModelDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const updateGroup = useGroupsStore((state) => state.updateGroup);

  if (!props.group) {
    return null;
  }

  const onSubmit = async (values: any) => {
    // form is valid
    const params: IGroupUpdate = {
      name: values.name,
      description: values.description,
      url: values.url,
      is_open: values.isOpen,
      tags: values.tags,
    };
    await updateGroup(props.group.id, params);
    props.handleClose();
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Edit Group"
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
              placeholder="Enter group name"
              defaultValue={props.group.name}
              inputRef={register({
                required: "Name is required",
              })}
            />
          </FormGroup>

          <FormGroup label="Description" labelFor="description-input">
            <InputGroup
              id="description-input"
              name="description"
              placeholder="Enter group description"
              defaultValue={props.group.description}
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
