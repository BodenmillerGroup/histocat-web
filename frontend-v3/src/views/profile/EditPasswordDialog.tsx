import { Button, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { IUserProfileUpdate } from "modules/profile/models";
import { useRef, useState } from "react";
import { useAuthStore } from "modules/auth";
import { Tooltip2 } from "@blueprintjs/popover2";

type EditUserDialogProps = {
  isOpen: boolean;
  handleClose(): void;
};

export function EditPasswordDialog(props: EditUserDialogProps) {
  const { register, errors, handleSubmit, watch } = useForm();
  const [showPassword, setShowPassword] = useState(false);
  const { updateUserPassword } = useAuthStore();

  const password = useRef({});
  password.current = watch("password", "");

  const onSubmit = async (values: IUserProfileUpdate) => {
    // form is valid
    await updateUserPassword(values.password!);
    props.handleClose();
  };

  const lockButton = (
    <Tooltip2 content={`${showPassword ? "Hide" : "Show"} Password`}>
      <Button
        icon={showPassword ? "unlock" : "lock"}
        intent={Intent.WARNING}
        minimal={true}
        onClick={() => setShowPassword(!showPassword)}
      />
    </Tooltip2>
  );

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Change Password"
      usePortal={true}
      isOpen={props.isOpen}
      className="bp3-dark"
      canOutsideClickClose={false}
    >
      <form onSubmit={handleSubmit(onSubmit)}>
        <div className={Classes.DIALOG_BODY}>
          <FormGroup
            label="Password"
            labelFor="password-input"
            labelInfo="(required)"
            intent="danger"
            helperText={errors?.password?.message}
          >
            <InputGroup
              id="password-input"
              name="password"
              placeholder="Enter password"
              type={showPassword ? "text" : "password"}
              rightElement={lockButton}
              inputRef={register({
                required: "Password is required",
                minLength: {
                  value: 3,
                  message: "Password must have at least 3 characters",
                },
              })}
            />
          </FormGroup>

          <FormGroup
            label="Confirm password"
            labelFor="password-confirm-input"
            labelInfo="(required)"
            intent="danger"
            helperText={errors?.passwordConfirm?.message}
          >
            <InputGroup
              id="password-confirm-input"
              name="passwordConfirm"
              placeholder="Confirm password"
              type={showPassword ? "text" : "password"}
              rightElement={lockButton}
              inputRef={register({
                required: "Password is required",
                validate: (value) => value === password.current || "The passwords do not match",
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
