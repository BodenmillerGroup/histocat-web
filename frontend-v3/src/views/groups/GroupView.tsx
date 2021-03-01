import { Route, Switch, useParams, useRouteMatch } from "react-router-dom";
import React, { useEffect, useState } from "react";
import { useProjectsStore } from "modules/projects";
import { LoadingView } from "components/LoadingView";
import { ProjectsView } from "./projects/ProjectsView";
import shallow from "zustand/shallow";
import { useGroupsStore } from "modules/groups";
import { ProjectContainer } from "./projects/ProjectContainer";

export default function GroupView() {
  const { path } = useRouteMatch();
  const { groupId } = useParams<{ groupId: string }>();
  const setActiveGroupId = useGroupsStore((state) => state.setActiveGroupId);
  const { getGroupProjects, getProjectsTags } = useProjectsStore(
    (state) => ({
      getGroupProjects: state.getGroupProjects,
      getProjectsTags: state.getProjectsTags,
    }),
    shallow
  );
  const [loading, setLoading] = useState<boolean>(true);

  // Mounted
  useEffect(() => {
    setActiveGroupId(+groupId);
    Promise.all([getGroupProjects(+groupId), getProjectsTags(+groupId)]).then(() => setLoading(false));
  }, [groupId, getGroupProjects, getProjectsTags, setActiveGroupId]);

  // Unmounted
  useEffect(() => {
    return () => {
      setActiveGroupId(null);
    };
  }, [setActiveGroupId]);

  if (loading) {
    return <LoadingView />;
  }

  return (
    <Switch>
      <Route exact={true} path={path} component={ProjectsView} />
      <Route path={`${path}/projects/:projectId`} component={ProjectContainer} />
    </Switch>
  );
}
