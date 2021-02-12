import { Button, Card, Elevation } from "@blueprintjs/core";
import { useProfileStore } from "../../modules/profile";
import styles from "./Profile.module.scss";
import { useHistory } from "react-router-dom";

export function ProfileView() {
  const userProfile = useProfileStore((state) => state.userProfile);
  const history = useHistory();

  return (
    <div className={styles.container}>
      <Card interactive={false} elevation={Elevation.TWO}>
        <h2>User Profile</h2>
        <h3>Name:</h3>
        <p>{userProfile?.name}</p>
        <h3>Email:</h3>
        <p>{userProfile?.email}</p>
        <div className={styles.actions}>
          <Button text="Edit profile" onClick={() => history.push("/profile/edit")} />
          <Button text="Change password" onClick={() => history.push("/profile/password")} />
        </div>
      </Card>
    </div>
  );
}
