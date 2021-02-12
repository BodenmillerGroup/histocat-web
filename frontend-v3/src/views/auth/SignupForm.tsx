import { useForm } from "react-hook-form";
import { Button, Card, Elevation, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useAuthStore } from "modules/auth";
import styles from "./Auth.module.scss";
import classNames from "classnames";
import { useRef, useState } from "react";
import { Tooltip2 } from "@blueprintjs/popover2";
import { appName } from "../../env";
import { AppToaster } from "utils/toaster";

export function SignupForm() {
  const { register, errors, handleSubmit, watch } = useForm();
  const { signUp, checkUserExists } = useAuthStore();
  const [showPassword, setShowPassword] = useState(false);

  const password = useRef({});
  password.current = watch("password", "");

  const onSubmit = async (values: any) => {
    // form is valid
    const userExist = await checkUserExists(values.email);
    if (userExist) {
      AppToaster.show({ message: "User with this email already exists", intent: "warning" });
    } else {
      await signUp({
        email: values.email,
        name: values.name,
        password: values.password,
      });
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
        <h3>{appName} - Signup</h3>
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

          <FormGroup label="Name" labelFor="name-input">
            <InputGroup id="name-input" name="name" placeholder="Enter your real name" inputRef={register({})} />
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
        </div>
      </Card>
    </div>
  );
}
