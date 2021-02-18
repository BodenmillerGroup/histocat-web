import { Button, Checkbox, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { IUserProfile, IUserProfileUpdate } from "modules/profile/models";
import { useUsersStore } from "modules/users";

type EditUserDialogProps = {
  isOpen: boolean;
  handleClose(): void;
  user: IUserProfile;
};

export function EditUserDialog(props: EditUserDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const updateUser = useUsersStore((state) => state.updateUser);

  const onSubmit = async (values: any) => {
    // form is valid
    const params: IUserProfileUpdate = {
      email: values.email,
      name: values.name,
      is_admin: values.isAdmin,
      is_active: values.isActive,
    };
    await updateUser(props.user.id, params);
    props.handleClose();
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Edit User"
      usePortal={true}
      isOpen={props.isOpen}
      className={Classes.DARK}
      canOutsideClickClose={false}
    >
      <form onSubmit={handleSubmit(onSubmit)}>
        <div className={Classes.DIALOG_BODY}>
          <FormGroup
            label="Email"
            labelFor="email-input"
            labelInfo="(required)"
            intent="danger"
            helperText={errors?.email?.message}
          >
            <InputGroup
              id="email-input"
              name="email"
              placeholder="Enter email"
              defaultValue={props.user.email}
              inputRef={register({
                required: "Email is required",
                pattern: {
                  value: /^[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}$/i,
                  message: "Invalid email address format",
                },
              })}
            />
          </FormGroup>

          <FormGroup label="Name" labelFor="name-input">
            <InputGroup
              id="name-input"
              name="name"
              placeholder="Enter your real name"
              defaultValue={props.user.name}
              inputRef={register({})}
            />
          </FormGroup>

          <Checkbox
            name="isActive"
            label="Active"
            defaultChecked={props.user.is_active}
            inline={true}
            inputRef={register({})}
          />
          <Checkbox
            name="isAdmin"
            label="Admin"
            defaultChecked={props.user.is_admin}
            inline={true}
            inputRef={register({})}
          />
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
