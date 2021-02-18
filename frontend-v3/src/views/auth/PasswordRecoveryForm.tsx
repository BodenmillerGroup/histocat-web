import { useForm } from "react-hook-form";
import { Button, Card, Classes, Elevation, FormGroup, InputGroup } from "@blueprintjs/core";
import { useAuthStore } from "modules/auth";
import styles from "./Auth.module.scss";
import classNames from "classnames";
import { appName } from "../../env";

export function PasswordRecoveryForm() {
  const { register, errors, handleSubmit } = useForm();
  const passwordRecovery = useAuthStore(state => state.passwordRecovery);

  const onSubmit = async (values: any) => {
    // form is valid
    await passwordRecovery(values.email);
  };

  return (
    <div className={styles.formContainer}>
      <Card interactive={false} elevation={Elevation.TWO} className={classNames(styles.form, Classes.DARK)}>
        <h3>{appName} - Password Recovery</h3>
        <h5>A password recovery email will be sent to the registered account</h5>
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
