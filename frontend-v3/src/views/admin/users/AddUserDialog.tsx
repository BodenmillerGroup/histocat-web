import { Button, Checkbox, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { IUserProfileCreate } from "modules/profile/models";
import { useUsersStore } from "modules/users";

type AddUserDialogProps = {
  isOpen: boolean;
  handleClose(): void;
};

export function AddUserDialog(props: AddUserDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const createUser = useUsersStore((state) => state.createUser);

  const onSubmit = async (values: any) => {
    // form is valid
    const params: IUserProfileCreate = {
      email: values.email,
      name: values.name,
      is_admin: values.isAdmin,
      is_active: values.isActive,
    };
    await createUser(params);
    props.handleClose();
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Add User"
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
            <InputGroup id="name-input" name="name" placeholder="Enter your real name" inputRef={register({})} />
          </FormGroup>

          <Checkbox name="isActive" label="Active" inline={true} defaultChecked={true} inputRef={register({})} />
          <Checkbox name="isAdmin" label="Admin" inline={true} inputRef={register({})} />
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
