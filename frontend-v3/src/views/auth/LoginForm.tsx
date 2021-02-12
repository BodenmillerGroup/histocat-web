import { useForm } from "react-hook-form";
import { Button, Callout, Card, Elevation, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useAuthStore } from "modules/auth";
import styles from "./Auth.module.scss";
import classNames from "classnames";
import { useState } from "react";
import { Tooltip2 } from "@blueprintjs/popover2";
import { appName } from "../../env";

export function LoginForm() {
  const { register, errors, handleSubmit } = useForm();
  const { login, loginError } = useAuthStore();
  const [showPassword, setShowPassword] = useState(false);

  const onSubmit = async (values: any) => {
    // form is valid
    await login(values.email, values.password);
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
        <h3>{appName} - Login</h3>
        <form onSubmit={handleSubmit(onSubmit)}>
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

          <Button type="submit" text="Submit" />
        </form>
        {loginError && (
          <Callout icon="warning-sign" intent="warning" className={styles.callout}>
            Incorrect email or password
          </Callout>
        )}
        <div className={styles.linkContainer}>
          <a href="/password-recovery" className={styles.link}>
            Recover password
          </a>
          <a href="/signup" className={styles.link}>
            Register an account
          </a>
        </div>
      </Card>
    </div>
  );
}
