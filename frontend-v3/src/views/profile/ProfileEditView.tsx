import { Button, Card, Elevation, FormGroup, InputGroup } from "@blueprintjs/core";
import { useProfileStore } from "../../modules/profile";
import styles from "./Profile.module.scss";
import { useForm } from "react-hook-form";
import { IUserProfileUpdate } from "modules/profile/models";
import { useHistory } from "react-router-dom";

export function ProfileEditView() {
  const { register, errors, handleSubmit } = useForm();
  const history = useHistory();
  const { userProfile, updateUserProfile } = useProfileStore();

  const onSubmit = async (values: IUserProfileUpdate) => {
    // form is valid
    await updateUserProfile({
      email: values.email,
      name: values.name,
    });
    history.goBack();
  };

  return (
    <div className={styles.container}>
      <Card interactive={false} elevation={Elevation.TWO}>
        <h2>Edit Profile</h2>
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
              defaultValue={userProfile?.email}
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
              defaultValue={userProfile?.name}
              inputRef={register({})}
            />
          </FormGroup>

          <Button type="reset" text="Reset" />
          <Button type="submit" text="Submit" intent="primary" style={{ marginLeft: "10px" }} />
        </form>
      </Card>
    </div>
  );
}
