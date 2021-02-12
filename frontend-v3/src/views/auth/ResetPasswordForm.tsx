import { useForm } from "react-hook-form";
import { Button, Card, Elevation, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useAuthStore } from "modules/auth";
import styles from "./Auth.module.scss";
import classNames from "classnames";
import { useRef, useState } from "react";
import { Tooltip2 } from "@blueprintjs/popover2";
import { appName } from "../../env";
import history from "utils/history";
import { AppToaster } from "utils/toaster";

export function ResetPasswordForm() {
  const { register, errors, handleSubmit, watch } = useForm();
  const { resetPassword } = useAuthStore();
  const [showPassword, setShowPassword] = useState(false);

  const password = useRef({});
  password.current = watch("password", "");

  const checkToken = () => {
    const token = history.location.search.replace("?token=", "");
    if (!token || token.length < 16) {
      AppToaster.show({ message: "No valid token provided in the URL", intent: "danger" });
      history.push("/password-recovery");
    } else {
      return token;
    }
  };

  const onSubmit = async (values: any) => {
    // form is valid
    const token = checkToken();
    if (token) {
      await resetPassword(values.password, token);
    }
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
    <div className={styles.formContainer}>
      <Card interactive={false} elevation={Elevation.TWO} className={classNames(styles.form, "bp3-dark")}>
        <h3>{appName} - Reset Password</h3>
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

          <Button type="submit" text="Submit" />
        </form>
        <div className={styles.linkContainer}>
          <a href="/login" className={styles.link}>
            Already have an account?
          </a>
          <a href="/signup" className={styles.link}>
            Register an account
          </a>
        </div>
      </Card>
    </div>
  );
}
