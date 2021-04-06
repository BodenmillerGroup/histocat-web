import React from "react";
import TitleInfo from "vitessce/components/TitleInfo";
import { BlendView } from "./BlendView";
import { CoordinationScopes } from "vitessce/types";

type BlendViewSubscriberProps = {
  coordinationScopes: CoordinationScopes;
  removeGridComponent(event: any): void;
};

export function BlendViewSubscriber(props: BlendViewSubscriberProps) {
  const { removeGridComponent } = props;
  return (
    <TitleInfo
      title="Image"
      removeGridComponent={removeGridComponent}
      isSpatial={true}
      theme="dark"
      isReady={true}
    >
      <BlendView />
    </TitleInfo>
  );
}
