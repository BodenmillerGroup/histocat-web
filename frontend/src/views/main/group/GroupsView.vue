<template>
  <v-container fluid>
    <v-row class="mt-6 mx-6">
      <v-text-field
        v-model="search"
        append-icon="mdi-magnify"
        label="Search"
        single-line
        hide-details
        clearable
        solo
        class="mt-1"
      />
      <v-select
        v-model="selectedTags"
        :items="tags"
        chips
        deletable-chips
        clearable
        label="Tags"
        multiple
        solo
        class="mt-1 ml-4"
      />
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" dark fab color="primary lighten-1" to="/main/groups/create" class="ml-4">
            <v-icon>mdi-plus</v-icon>
          </v-btn>
        </template>
        <span>Create group</span>
      </v-tooltip>
    </v-row>
    <masonry :cols="{ default: 4, 1000: 3, 700: 2, 400: 1 }" :gutter="{ default: '0px' }">
      <GroupCard v-for="group in groups" :key="group.id" :group="group" :user="userProfile" />
    </masonry>
  </v-container>
</template>

<script lang="ts">
import GroupCard from "@/views/main/group/GroupCard.vue";
import { mainModule } from "@/modules/main";
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";

@Component({
  components: { GroupCard },
})
export default class GroupsView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly groupContext = groupModule.context(this.$store);

  search = "";
  selectedTags: string[] = [];

  get tags(): any[] {
    const list = this.groupContext.getters.tags;
    return list.map((item) => {
      return {
        text: item,
      };
    });
  }

  get userProfile() {
    return this.mainContext.getters.userProfile;
  }

  get groups() {
    const items = this.search
      ? this.groupContext.getters.groups.filter((item) => item.name.toLowerCase().includes(this.search.toLowerCase()))
      : this.groupContext.getters.groups;
    if (this.selectedTags.length === 0) {
      return items;
    } else {
      return items.filter((group) => {
        if (group.tags && this.selectedTags.some((r) => group.tags.includes(r))) {
          return group;
        }
      });
    }
  }

  async mounted() {
    await Promise.all([this.groupContext.actions.getGroups(), this.groupContext.actions.getTags()]);
  }
}
</script>
