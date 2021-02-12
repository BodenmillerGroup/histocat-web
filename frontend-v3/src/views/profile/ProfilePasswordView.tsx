import { Button, Card, Elevation, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import styles from "./Profile.module.scss";
import { useForm } from "react-hook-form";
import { useHistory } from "react-router-dom";
import { useRef, useState } from "react";
import { Tooltip2 } from "@blueprintjs/popover2";
import { IUserProfileUpdate } from "../../modules/profile/models";
import { useAuthStore } from "../../modules/auth";

export function ProfilePasswordView() {
  const { register, errors, handleSubmit, watch } = useForm();
  const history = useHistory();
  const [showPassword, setShowPassword] = useState(false);
  const { updateUserPassword } = useAuthStore();

  const password = useRef({});
  password.current = watch("password", "");

  const onSubmit = async (values: IUserProfileUpdate) => {
    // form is valid
    await updateUserPassword(values.password!);
    history.goBack();
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
    <div className={styles.container}>
      <Card interactive={false} elevation={Elevation.TWO}>
        <h2>Change Password</h2>
        <form onSubmit={handleSubmit(onSubmit)}>
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

          <Button type="reset" text="Reset" />
          <Button type="submit" text="Submit" intent="primary" style={{ marginLeft: "10px" }} />
        </form>
      </Card>
    </div>
  );
}
