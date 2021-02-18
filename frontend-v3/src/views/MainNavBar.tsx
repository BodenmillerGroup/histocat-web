import {
  Alignment,
  Button,
  Menu,
  MenuItem,
  Navbar,
  NavbarDivider,
  NavbarGroup,
  NavbarHeading,
} from "@blueprintjs/core";
import { Popover2 } from "@blueprintjs/popover2";
import { useAuthStore } from "modules/auth";
import { appName } from "../env";
import { useHistory } from "react-router-dom";
import { useGroupsStore } from "modules/groups";
import { useProfileStore } from "modules/profile";
import { LayoutsControl } from "./LayoutsControl";

export function MainNavBar() {
  const hasAdminAccess = useProfileStore((state) => state.hasAdminAccess());
  const userLogout = useAuthStore((state) => state.userLogout);
  const activeGroupId = useGroupsStore((state) => state.activeGroupId);
  const history = useHistory();

  return (
    <Navbar fixedToTop={false}>
      <NavbarGroup align={Alignment.LEFT}>
        <NavbarHeading>
          <Button minimal={true} text={appName} onClick={() => history.push("/")} />
        </NavbarHeading>
        <NavbarDivider />
        {activeGroupId && (
          <Button minimal={true} icon="projects" text="Projects" onClick={() => history.push(`/${activeGroupId}`)} />
        )}
        <LayoutsControl />
      </NavbarGroup>
      <NavbarGroup align={Alignment.RIGHT}>
        {hasAdminAccess && (
          <Popover2
            content={
              <Menu>
                <MenuItem text="Users" icon="user" onClick={() => history.push("/admin/users")} />
                <MenuItem text="Models" icon="predictive-analysis" onClick={() => history.push("/admin/models")} />
              </Menu>
            }
            placement="bottom"
          >
            <Button minimal={true} icon="key" text="Admin" />
          </Popover2>
        )}
        <Popover2
          content={
            <Menu>
              <MenuItem text="Profile" href="/profile" />
              <MenuItem text="Logout" onClick={userLogout} />
            </Menu>
          }
          placement="bottom"
        >
          <Button minimal={true} icon="user" />
        </Popover2>
        <Button minimal={true} icon="notifications" />
        <Button minimal={true} icon="cog" />
      </NavbarGroup>
    </Navbar>
  );
}
