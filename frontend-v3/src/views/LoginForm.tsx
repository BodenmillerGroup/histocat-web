import { useForm } from "react-hook-form";
import { Button, Card, Elevation, FormGroup, InputGroup } from "@blueprintjs/core";
import { useAuthStore } from "state/auth";
import { useHistory } from "react-router";

export function LoginForm() {
  const { register, errors, handleSubmit } = useForm();
  const { login } = useAuthStore();

  const onSubmit = async (values: any) => {
    // form is valid
    login(values.email, values.password);
    // eslint-disable-next-line react-hooks/rules-of-hooks
    const history = useHistory();
    history.push("/main")
  };

  return (
    <Card interactive={false} elevation={Elevation.TWO} className="bp3-dark">
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
            inputRef={register({
              required: "Password is required",
              validate: (value) => value.length > 2 || "Password must be 3 characters at minimum",
            })}
          />
        </FormGroup>

        <Button type="submit">Submit</Button>
      </form>
    </Card>
  );
}
