import { Button, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { IUserProfileUpdate } from "modules/profile/models";
import { useProfileStore } from "modules/profile";

type EditUserDialogProps = {
  isOpen: boolean;
  handleClose(): void;
};

export function EditProfileDialog(props: EditUserDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const userProfile = useProfileStore((state) => state.userProfile);
  const updateUserProfile = useProfileStore((state) => state.updateUserProfile);

  if (!userProfile) {
    return null;
  }

  const onSubmit = async (values: any) => {
    // form is valid
    const params: IUserProfileUpdate = {
      email: values.email,
      name: values.name,
    };
    await updateUserProfile(params);
    props.handleClose();
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Edit Profile"
      usePortal={true}
      isOpen={props.isOpen}
      className="bp3-dark"
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
              defaultValue={userProfile!.email}
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
              defaultValue={userProfile!.name}
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
